% version 202401

%% Define the number of points in the colormap
n = 256; % You can adjust this number based on your preference

% Create a custom colormap
% Linearly interpolate between blue (at the start), white (in the middle), and red (at the end)
blueToWhite = [linspace(0, 1, n/2)' linspace(0, 1, n/2)' ones(n/2, 1)];
whiteToRed = [ones(n/2, 1) linspace(1, 0, n/2)' linspace(1, 0, n/2)'];
customColormap = [blueToWhite; whiteToRed(2:end,:)];

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

%       load('Supriya_data_metrics.mat')
%         supriya_index = [4;27;22;6;15];
%         supriya_dual_ligand_index = [3;5;8;9;26;14;25;18;19];

%         % sampled data
%         data_filename = 'Sim5_codon_all5dose_metric.mat';
%         load(strcat(data_save_file_path_1,data_filename))
%         sample_dual_ligand_index = [1;2;3;4;10;9;6;8;5;7];
%         sample_index = [13;12;11;14;15];


%% draw multi ligand stim match
if 1 % draw all comb ligand stim for weighted
    fig_save_path = '../SubFigures2023/';
    
    cal_responder_ratio = 1;
    integral_pos_threshold_vec = [0.2,0.3,0.4,0.5,0.6,0.7,0.8];
    peak_threshold_vec = [0.08,0.10,0.12,0.14,0.16,0.18,0.20];
    
    % Use arrayfun for element-wise operation
    integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
    peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);
    
    if cal_responder_ratio
        figure(1)
        paperpos = [0,0,220,50]*1.8;
        papersize = [220,50]*1.8;
        draw_pos = [10,10,200,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
         
        for i_thresh = 1:length(integral_pos_threshold_vec)
            integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
            peak_threshold = peak_threshold_vec(i_thresh);
            
            load('../raw_data2023/Sim10_all_single_Comb_ligands_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            % 0.25 * 2
            % 0.1 * 5
            
            load('../raw_data2023/Sim9_all_Comb_ligands_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(5+i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(5+i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            responder_ratio_sim{i_thresh} = responder_ratio;
            responder_ratio_by_peak_sim{i_thresh} = responder_ratio_by_peak;
            
%             figure(1)
%             plot(1:length(responder_ratio),responder_ratio,'LineWidth',2);hold on
            
             
        end
        
%         Xtick_labels ={};
%         for i_Xtick_labels = 1:length(responder_ratio)
%             Xtick_labels(i_Xtick_labels) = {''};
%         end
%         
%         figure(1)
%         %legend(integral_pos_threshold_legend)
%         
%         set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
%             'FontSize',8)
%         ylim([0,1])
%         xlim([0.5,31.5])
        %saveas(gcf,strcat(fig_save_path,'Responder_ratio_by_integral_sim_weighted_match'),'epsc')
        %close()
        
    end
    
    %%
    analysis_stad = 0;
    if analysis_stad
        figure(1)
        hist(data.model_sim{1, 1}(1:9:9000,end),20)
        title('end point')
        xlabel('nNFkB')
        ylabel('counts')
        
        figure(2)
        hist(data.model_sim{1, 1}(1:9:9000,1),20)
        title('start point')
        xlabel('nNFkB')
        ylabel('counts')
    end
    
    %%
    plot_heatmap = 0;
    if plot_heatmap
        
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 2997;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 integrals_pos 97
        metric_order_name_{1} = 'integrals_pos';
        metric_order_column = 97;
        
        order_name_vec = {'ligand3_stim'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1});
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        for i_cond = 1:length(data_NFkB.info_ligand)
            [~,data_NFkB.order{i_cond}] = sort(metrics{i_cond}.(metric_order_name_{1})(1:9:9000,metric_order_column),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{1})
        
    end
end

%% draw multi ligand stim weighted
if 0 % draw all comb ligand stim for weighted
    fig_save_path = '../SubFigures2023/';
    
    cal_responder_ratio = 1;
    integral_pos_threshold_vec = [0.2,0.3,0.4,0.5,0.6,0.7,0.8];
    peak_threshold_vec = [0.08,0.10,0.12,0.14,0.16,0.18,0.20];
    
    % Use arrayfun for element-wise operation
    integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
    peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);
    
    if cal_responder_ratio
        figure(1)
        paperpos = [0,0,220,50]*1.8;
        papersize = [220,50]*1.8;
        draw_pos = [10,10,200,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
         
        for i_thresh = 1:length(integral_pos_threshold_vec)
            integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
            peak_threshold = peak_threshold_vec(i_thresh);
            
            load('../raw_data2023/Sim10_all_single_Comb_ligands_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            % 0.25 * 2
            % 0.1 * 5
            
            load('../raw_data2023/Sim9_all_Comb_ligands_weighted_match_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(5+i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(5+i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            responder_ratio_sim{i_thresh} = responder_ratio;
            responder_ratio_by_peak_sim{i_thresh} = responder_ratio_by_peak;
            
            figure(1)
            plot(1:length(responder_ratio),responder_ratio,'LineWidth',2);hold on
            
             
        end
        
        Xtick_labels ={};
        for i_Xtick_labels = 1:length(responder_ratio)
            Xtick_labels(i_Xtick_labels) = {''};
        end
        
        figure(1)
        %legend(integral_pos_threshold_legend)
        
        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)
        ylim([0,1])
        xlim([0.5,31.5])
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_by_integral_sim_weighted_match'),'epsc')
        close()
        
    end
    
    %%
    analysis_stad = 0;
    if analysis_stad
        figure(1)
        hist(data.model_sim{1, 1}(1:9:9000,end),20)
        title('end point')
        xlabel('nNFkB')
        ylabel('counts')
        
        figure(2)
        hist(data.model_sim{1, 1}(1:9:9000,1),20)
        title('start point')
        xlabel('nNFkB')
        ylabel('counts')
    end
    
    %%
    plot_heatmap = 0;
    if plot_heatmap
        
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 2997;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 integrals_pos 97
        metric_order_name_{1} = 'integrals_pos';
        metric_order_column = 97;
        
        order_name_vec = {'ligand3_stim'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1});
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        for i_cond = 1:length(data_NFkB.info_ligand)
            [~,data_NFkB.order{i_cond}] = sort(metrics{i_cond}.(metric_order_name_{1})(1:9:9000,metric_order_column),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{1})
        
    end
end

%% draw multi ligand stim for scrambling parameters
if 0 % draw all comb ligand stim
    fig_save_path = '../SubFigures2023/';
    
    cal_responder_ratio = 1;
    integral_pos_threshold_vec = [0.2,0.3,0.4,0.5,0.6,0.7,0.8];
    peak_threshold_vec = [0.08,0.10,0.12,0.14,0.16,0.18,0.20];
    
    % Use arrayfun for element-wise operation
    integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
    peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);
    
    if cal_responder_ratio
        figure(1)
        paperpos = [0,0,220,50]*1.8;
        papersize = [220,50]*1.8;
        draw_pos = [10,10,200,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        figure(2)
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        
        for i_thresh = 1:length(integral_pos_threshold_vec)
            integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
            peak_threshold = peak_threshold_vec(i_thresh);
            
            load('../raw_data2023/Sim10_all_single_Comb_ligands_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            % 0.25 * 2
            % 0.1 * 5
            
            load('../raw_data2023/Sim11_all_comb_scramble_codon_metric.mat')
            
            for i_cond = 1:length(metrics)
                responder_ratio(5+i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:end,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:end,end));
                responder_ratio_by_peak(5+i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:end,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:end,end));
                
            end
            responder_ratio_sim{i_thresh} = responder_ratio;
            responder_ratio_by_peak_sim{i_thresh} = responder_ratio_by_peak;
            
            figure(1)
            plot(1:length(responder_ratio),responder_ratio,'LineWidth',2);hold on
            
            figure(2)
            plot(1:length(responder_ratio_by_peak),responder_ratio_by_peak,'LineWidth',2);hold on
            
        end
        
        Xtick_labels ={};
        for i_Xtick_labels = 1:length(responder_ratio)
            Xtick_labels(i_Xtick_labels) = {''};
        end
        
        figure(1)
        %legend(integral_pos_threshold_legend)
        
        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)
        ylim([0,1])
        xlim([0.5,31.5])
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_by_integral_sim_scramble'),'epsc')
        close()
        
        figure(2)
        %legend(peak_threshold_legend)
        set(gca,'XTick',1:length(responder_ratio),'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
            'FontSize',8)
        
        ylim([0,1])
        xlim([0.5,31.5])
        
        saveas(gcf,strcat(fig_save_path,'Responder_ratio_by_peak_sim_scramble'),'epsc')
        close()
    end
    
    %%
    analysis_stad = 0;
    if analysis_stad
        figure(1)
        hist(data.model_sim{1, 1}(1:9:9000,end),20)
        title('end point')
        xlabel('nNFkB')
        ylabel('counts')
        
        figure(2)
        hist(data.model_sim{1, 1}(1:9:9000,1),20)
        title('start point')
        xlabel('nNFkB')
        ylabel('counts')
    end
    
    %%
    plot_heatmap = 0;
    if plot_heatmap
        
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 2997;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 integrals_pos 97
        metric_order_name_{1} = 'integrals_pos';
        metric_order_column = 97;
        
        order_name_vec = {'ligand3_stim'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1});
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        for i_cond = 1:length(data_NFkB.info_ligand)
            [~,data_NFkB.order{i_cond}] = sort(metrics{i_cond}.(metric_order_name_{1})(1:9:9000,metric_order_column),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{1})
        
    end
end

%% compare with exp. multi ligand
if 1 % compare with exp. 3 ligand
    load('Supriya_data_metrics_2024.mat')
    
    cal_responder_ratio = 1;
    index_exp_ligand_vects = {[4];% 1-TNF
        [6,24]; %2-
        [1,10]; %3
        [7,12,15,17,22]; %4
        [2,11,23,27]; %5
        [8];%6
        [5];%7
        [9];%8
        [3];%9
        [25];%10
        [19];%11
        [26];%12
        [13,16,18];%13
        [];%14
        [14];%15
        [];%16
        [];%17
        [];%18
        [20];%19
        [];%20
        [21];%21
        [];%22
        [];%23
        [];%24
        [];%25
        [];%26
        [];%27
        [];%28
        [];%29
        [];%30
        [28]};%31
    if cal_responder_ratio
        
        %% scatter plot
        scatter_plot = 1;
        if scatter_plot
            integral_pos_threshold = 0.3;
            peak_threshold = 0.15;
            % 0.25 * 2
            % 0.1 * 5
            
            clear responder_ratio
            index_scatter_y = 1;
            for i_cond = 1:length(index_exp_ligand_vects)
                index_exp_ligand = index_exp_ligand_vects{i_cond};
                responder_ratio(i_cond) = 0;
                responder_ratio_by_peak(i_cond) = 0;
                for i_rpc = 1:length(index_exp_ligand)
                    responder_ratio_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end) > integral_pos_threshold)/length(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end));
                    responder_ratio_scatter_x(index_scatter_y) = i_cond;
                    responder_ratio(i_cond) =responder_ratio(i_cond)+ responder_ratio_scatter_y(index_scatter_y);
                    
                    responder_ratio_by_peak_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end) > peak_threshold)/length(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end));
                    responder_ratio_by_peak_scatter_x(index_scatter_y) = i_cond;
                    responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)+responder_ratio_by_peak_scatter_y(index_scatter_y);
                    index_scatter_y = index_scatter_y+1;
                end
                
                if responder_ratio(i_cond)
                    responder_ratio(i_cond) = responder_ratio(i_cond)/length(index_exp_ligand);
                    responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)/length(index_exp_ligand);
                end
            end
            
            figure(1)
            responder_ratio_nonzero_index = find(responder_ratio > 0);
            plot(responder_ratio_nonzero_index,responder_ratio(responder_ratio_nonzero_index));hold on
            scatter(responder_ratio_scatter_x,responder_ratio_scatter_y,'filled');hold
            
            figure(2)
            responder_ratio_by_peak_nonzero_index = find(responder_ratio > 0);
            plot(responder_ratio_by_peak_nonzero_index,responder_ratio_by_peak(responder_ratio_by_peak_nonzero_index));hold on
            scatter(responder_ratio_by_peak_scatter_x,responder_ratio_by_peak_scatter_y,'filled');hold
            
        end
        
        %% plot lines
        plot_lines = 1;
        if plot_lines
            figure(1)
            paperpos = [0,0,220,50]*1.8;
            papersize = [220,50]*1.8;
            draw_pos = [10,10,200,30]*1.8;
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            figure(2)
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            
            
            integral_pos_threshold_vec = [0.05,0.1,0.15,0.2,0.25,0.3,0.40];
            peak_threshold_vec = [0.08,0.10,0.12,0.14,0.16,0.18,0.20];
            
            % Use arrayfun for element-wise operation
            integral_pos_threshold_legend = arrayfun(@num2str, integral_pos_threshold_vec, 'UniformOutput', false);
            peak_threshold_legend = arrayfun(@num2str, peak_threshold_vec, 'UniformOutput', false);
            
            for i_thresh = 1:length(integral_pos_threshold_vec)
                integral_pos_threshold = integral_pos_threshold_vec(i_thresh);
                peak_threshold = peak_threshold_vec(i_thresh);
                
                clear responder_ratio responder_ratio_scatter_y responder_ratio_scatter_x
                clear responder_ratio_by_peak_scatter_y responder_ratio_by_peak_scatter_x responder_ratio_by_peak
                index_scatter_y = 1;
                for i_cond = 1:length(index_exp_ligand_vects)
                    index_exp_ligand = index_exp_ligand_vects{i_cond};
                    responder_ratio(i_cond) = 0;
                    responder_ratio_by_peak(i_cond) = 0;
                    for i_rpc = 1:length(index_exp_ligand)
                        responder_ratio_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end) > integral_pos_threshold)/length(metrics{index_exp_ligand(i_rpc)}.integrals_pos(:,end));
                        responder_ratio_scatter_x(index_scatter_y) = i_cond;
                        responder_ratio(i_cond) =responder_ratio(i_cond)+ responder_ratio_scatter_y(index_scatter_y);
                        
                        responder_ratio_by_peak_scatter_y(index_scatter_y) = sum(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end) > peak_threshold)/length(metrics{index_exp_ligand(i_rpc)}.max_amplitude(:,end));
                        responder_ratio_by_peak_scatter_x(index_scatter_y) = i_cond;
                        responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)+responder_ratio_by_peak_scatter_y(index_scatter_y);
                        index_scatter_y = index_scatter_y+1;
                    end
                    
                    if responder_ratio(i_cond)
                        responder_ratio(i_cond) = responder_ratio(i_cond)/length(index_exp_ligand);
                        responder_ratio_by_peak(i_cond) = responder_ratio_by_peak(i_cond)/length(index_exp_ligand);
                    end
                end
                
                responder_ratio_exp{i_thresh} = responder_ratio;
                responder_ratio_scatter_y_exp{i_thresh} = responder_ratio_scatter_y;
                responder_ratio_scatter_x_exp{i_thresh} = responder_ratio_scatter_x;
                responder_ratio_by_peak_exp{i_thresh} = responder_ratio_by_peak;
                
                figure(1)
                responder_ratio_nonzero_index = find(responder_ratio > 0);
                plot(responder_ratio_nonzero_index,responder_ratio(responder_ratio_nonzero_index),'LineWidth',2);hold on
                % scatter(responder_ratio_scatter_x,responder_ratio_scatter_y,'filled');hold
                
            end
            
            length_responder_ratio_total = 31;
            Xtick_labels ={};
            for i_Xtick_labels = 1:length_responder_ratio_total
                Xtick_labels(i_Xtick_labels) = {''};
            end
            
            figure(1)
            %legend(integral_pos_threshold_legend)
            
            set(gca,'XTick',1:length_responder_ratio_total,'XTickLabel',Xtick_labels,'YTick',[0:0.2:1],'YTickLabel',{},...
                'FontSize',8)
            ylim([0,1])
            xlim([0.5,31.5])
            saveas(gcf,strcat(fig_save_path,'Responder_ratio_by_integral_exp_with7sti'),'epsc')
            close()

            
        end
        
        %% compare exp with sim
        compare_exp_sim = 1;
        if compare_exp_sim
            to_put_into_draw_3ligand_compare_multi_stim_exp_sim
        end
        
    end
    
   
end

%% 5 ligand stim
if 0 % 5 ligand stim
    load('../raw_data2023/Sim7_5ligands_codon_metric.mat')
    
    cal_responder_ratio = 1;
    if cal_responder_ratio
        integral_pos_threshold = 0.3;
        peak_threshold = 0.15;
        % 0.25 * 2
        % 0.1 * 5
        for i_cond = 1:length(metrics)
            responder_ratio(i_cond) = sum(metrics{i_cond}.integrals_pos(1:9:9000,end) > integral_pos_threshold)/length(metrics{i_cond}.integrals_pos(1:9:9000,end));
            responder_ratio_by_peak(i_cond) = sum(metrics{i_cond}.max_amplitude(1:9:9000,end) > peak_threshold)/length(metrics{i_cond}.max_amplitude (1:9:9000,end));
        end
    end
    
    analysis_stad = 0;
    if analysis_stad
        figure(1)
        hist(data.model_sim{1, 1}(1:9:9000,end),20)
        title('end point')
        xlabel('nNFkB')
        ylabel('counts')
        
        figure(2)
        hist(data.model_sim{1, 1}(1:9:9000,1),20)
        title('start point')
        xlabel('nNFkB')
        ylabel('counts')
    end
    
    plot_heatmap = 0;
    if plot_heatmap
        
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 2997;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 integrals_pos 97
        metric_order_name_{1} = 'integrals_pos';
        metric_order_column = 97;
        
        order_name_vec = {'ligand5_stim'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1});
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        for i_cond = 1:length(data_NFkB.info_ligand)
            [~,data_NFkB.order{i_cond}] = sort(metrics{i_cond}.(metric_order_name_{1})(1:9:9000,metric_order_column),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{1})
        
    end
end

%% draw 5 ligand single ligand stim

if 0 %draw 5 ligand single ligand stim
    load('../raw_data2023/Sim8_5_signle_ligand_codon_metric.mat')
    
    %% draw codons
    if 1
        % codon_mat = double;
        kmeans_cluster_plots = 0;
        single_cell_specificity_plot = 1;
        
        for i_data = 1:length(data.model_sim)
            
            data_info.info_ligand{i_data} = data.info_ligand{i_data};
            data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
            data_info.data_label{i_data} = 'SamplingComb';
            
        end
        metrics_cal = cell(1,length(metrics));
        metric_names = fieldnames(metrics{1});
        
        for i_metric_index = 1:length(metrics)
            for i_metric_name = 1:length(metric_names)
                metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
            end
        end
        
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        % codon_list = {'TotalActivity'};
        
        
        
        if kmeans_cluster_plots
            
            clear codon_mat
            i_column = 1;
            for i_ligand = 1:length(collect_feature_vects.info_ligand)
                for i_codon = 1:length(codon_list)
                    codon_mat(:,i_column) = collect_feature_vects.(codon_list{i_codon}){i_ligand}(:,:);
                    i_column = i_column + 1;
                end
            end
            % 3 & 3
            ligand_num = 5;
            
            cluster_num = 5;
            ra_seed_num = 3;
            
            %[~,order_inx] = sort(codon_mat(:,1),'ascend');
            
            figure(1)
            xlabels1 = cell(1,size(codon_mat,2));
            i_xlabels = 1;
            for i_ligand = 1:ligand_num
                for i_codon = 1:length(codon_list)
                    xlabels1{i_xlabels} = strcat(collect_feature_vects.info_ligand{i_ligand}(1:3),'-',codon_list{i_codon}(1));
                    i_xlabels = i_xlabels +1;
                end
            end
            h1 = heatmap(xlabels1,1:size(codon_mat,1),codon_mat(:,:),'GridVisible','off','Colormap',customColormap);%order_inx
            caxis([-1.5,1.5])
            
            YLabels = 1:size(codon_mat,1);
            % Convert each number in the array into a string
            YCustomXLabels1 = string(YLabels);
            % Replace all but the fifth elements by spaces
            YCustomXLabels1(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h1.YDisplayLabels = YCustomXLabels1;
            
            
            %         rng(ra_seed_num)
            %         inx = kmeans(codon_mat(:,[2,5]),cluster_num);
            
            rng(ra_seed_num)
            inx = kmeans(codon_mat(:,:),cluster_num);
            [inx_ordered,order_inx] = sort(inx,'ascend');
            
            i_cluster =1;
            Num_i_cluster = 1:sum(inx_ordered==i_cluster);
            for i_cluster = 2:cluster_num
                Num_i_cluster = [Num_i_cluster,(Num_i_cluster(end)+2):(Num_i_cluster(end)+2+sum(inx_ordered==i_cluster)-1)];
            end
            
            Num_i_ligand = [1:6,8:13,15:20,22:27,29:34];
            codon_to_plot = NaN(size(codon_mat,1)+cluster_num-1,size(codon_mat,2)+ligand_num-1);
            
            codon_to_plot(Num_i_cluster,Num_i_ligand) = codon_mat(order_inx,:);
            
            ylabels = {};
            xlabels = cell(1,size(codon_to_plot,2));
            i_xlabels = 1;
            for i_ligand = 1:ligand_num
                for i_codon = 1:length(codon_list)
                    xlabels{i_xlabels} = strcat(collect_feature_vects.info_ligand{i_ligand}(1:3),'-',codon_list{i_codon}(1));
                    i_xlabels = i_xlabels +1;
                end
                xlabels{i_xlabels} = collect_feature_vects.info_ligand{i_ligand}(1:3);
                i_xlabels = i_xlabels +1;
            end
            
            YLabels = 1:size(codon_to_plot,1);
            % Convert each number in the array into a string
            YCustomXLabels = string(YLabels);
            % Replace all but the fifth elements by spaces
            %YCustomXLabels(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            %h.YDisplayLabels = YCustomXLabels;
            
            figure(2)
            
            h2 = heatmap(xlabels(1:end-1),1:size(codon_to_plot,1),codon_to_plot,'GridVisible','off','Colormap',customColormap);%
            caxis([-1.5,1.5])
            
            title({strcat('#cluster=',num2str(cluster_num),', rng seed=',num2str(ra_seed_num))});
            
            
            % Replace all but the fifth elements by spaces
            YCustomXLabels(:) = " ";
            % Set the 'XDisplayLabels' property of the heatmap
            % object 'h' to the custom x-axis tick labels
            h2.YDisplayLabels = YCustomXLabels;
        end
        
        if single_cell_specificity_plot
            clear codon_single_cell
            codon_single_cell=cell(2,1);
            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:length(collect_feature_vects.info_ligand)
                    for i_codon = 1:length(codon_list)
                        codon_single_cell{i_cell}(i_ligand,i_codon) = collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
                    end
                end
            end
            
            clusterResults = cell2mat(cellfun(@applyDBSCAN, codon_single_cell, 'UniformOutput', false));
            
            %confused traj
            % espilon = 1 => 668
            
            cell_idx = 58;
            figure(1)
            for i_ligand = 1:5
                plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
            end
            ylabel('nNFkB (\mu M)')
            xlabel('Time (5min)')
            ylim([0,0.4])
            legend(collect_feature_vects.info_ligand)
            set(gca,'fontweight','b')
            %close()
            %specific traj
            % espilon = 2 => 709
            
            % espilon = 3 => 709
            cell_idx = 752;
            figure(2)
            for i_ligand = 1:5
                plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
            end
            ylabel('nNFkB (\mu M)')
            xlabel('Time (5min)')
            legend(collect_feature_vects.info_ligand)
            set(gca,'fontweight','b')
            close()
        end
    end
    
    
    %% draw heat map
    if 0
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1
        metric_order_name_{1} = 'oscpower';
        metric_order_column = 1;
        
        order_name_vec = {'sc_srs_TNF','sc_srs_LPS','sc_srs_CpG','sc_srs_PolyIC','sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1});
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 4% 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{i_order_i_ligand}.(metric_order_name_{1})(1:9:end,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2023_05(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
end

%% document epsilon for SRS
if 0 %draw 5 ligand single ligand stim
    load('../raw_data2023/Sim8_5_signle_ligand_codon_metric.mat')
    
    %% draw codons
    if 1
        % codon_mat = double;
        kmeans_cluster_plots = 0;
        single_cell_specificity_plot = 1;
        
        for i_data = 1:length(data.model_sim)
            
            data_info.info_ligand{i_data} = data.info_ligand{i_data};
            data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
            data_info.data_label{i_data} = 'SamplingComb';
            
        end
        metrics_cal = cell(1,length(metrics));
        metric_names = fieldnames(metrics{1});
        
        for i_metric_index = 1:length(metrics)
            for i_metric_name = 1:length(metric_names)
                metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
            end
        end
        
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        % codon_list = {'TotalActivity'};
        
        
        if single_cell_specificity_plot
            clear codon_single_cell
            codon_single_cell=cell(2,1);
            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:length(collect_feature_vects.info_ligand)
                    for i_codon = 1:length(codon_list)
                        codon_single_cell{i_cell}(i_ligand,i_codon) = collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
                    end
                end
            end
            
            clusterResults = cell2mat(cellfun(@applyDBSCAN, codon_single_cell, 'UniformOutput', false));
            
            %confused traj
            % espilon = 1 => 668
            
            cell_idx_vec = [668,58,752,709];
            clusterResults(cell_idx_vec)
            if 0
                for i_cell_idx = 1:length(cell_idx_vec)
                    cell_idx = cell_idx_vec(i_cell_idx);
                    figure(1)
                    paperpos = [0,0,70,50]*1.8;
                    papersize = [70,50]*1.8;
                    draw_pos = [10,10,50,30]*1.8;
                    
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    
                    for i_ligand = 1:5
                        plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
                    end
                    ylabel('nNFkB (\mu M)')
                    xlabel('Time (5min)')
                    ylim([0,0.4])
                    % legend(collect_feature_vects.info_ligand)
                    set(gca,'fontweight','b')
                    saveas(gcf,strcat(fig_save_path,'SRS_traj_',num2str(cell_idx)),'epsc')
                    close()
                end
            end
            %specific traj
            % espilon = 2 => 709
            
            % espilon = 3 => 709
            %             cell_idx = 752;
            %             figure(2)
            %             for i_ligand = 1:5
            %                 plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
            %             end
            %             ylabel('nNFkB (\mu M)')
            %             xlabel('Time (5min)')
            %             legend(collect_feature_vects.info_ligand)
            %             set(gca,'fontweight','b')
            %             close()
            %         end
        end
        
        
    end
end

%% pairwise SRS, cluter, MI
if 0 %draw 5 ligand single ligand stim
    load('../raw_data2023/Sim8_5_signle_ligand_codon_metric.mat')
    
    %% draw codons
    if 1
        % codon_mat = double;
        kmeans_cluster_plots = 0;
        single_cell_specificity_plot = 1;
        
        for i_data = 1:length(data.model_sim)
            
            data_info.info_ligand{i_data} = data.info_ligand{i_data};
            data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
            data_info.data_label{i_data} = 'SamplingComb';
            
        end
        metrics_cal = cell(1,length(metrics));
        metric_names = fieldnames(metrics{1});
        
        for i_metric_index = 1:length(metrics)
            for i_metric_name = 1:length(metric_names)
                metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
            end
        end
        
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        % codon_list = {'TotalActivity'};
        
        
        if single_cell_specificity_plot
            clear codon_single_cell
            pairwise_ligand_index = {[1,2];[1,3];[1,4];[1,5];[2,3];[2,4];[2,5];[3,4];[3,5];[4,5]};
            codon_single_cell=cell(2,length(pairwise_ligand_index));
            
            for j_cell = 1:length(pairwise_ligand_index)
                for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                    for i_ligand = 1:length(pairwise_ligand_index{j_cell})
                        for i_codon = 1:length(codon_list)
                            codon_single_cell{i_cell,j_cell}(i_ligand,i_codon) = collect_feature_vects.(codon_list{i_codon}){pairwise_ligand_index{j_cell}(i_ligand)}(i_cell,:);
                        end
                    end
                end
            end
            
            clusterResults = cell2mat(cellfun(@applyDBSCAN, codon_single_cell, 'UniformOutput', false));
            
            MI_bars = mean(log2(clusterResults));
            bar(MI_bars)
        end
        
        
    end
end


%% pairwise SRS, cluter, MI, dual ligands
if 0 %draw 5 ligand single ligand stim''

    load('../raw_data2023/Sim12_dual_signle_ligand_codon_metric.mat')
    
    %% draw codons
    if 1
        % codon_mat = double;
        kmeans_cluster_plots = 0;
        single_cell_specificity_plot = 1;
        
        for i_data = 1:length(data.model_sim)
            
            data_info.info_ligand{i_data} = data.info_ligand{i_data};
            data_info.info_dose_str{i_data} = data.info_dose_str{i_data};
            data_info.data_label{i_data} = 'SamplingComb';
            
        end
        metrics_cal = cell(1,length(metrics));
        metric_names = fieldnames(metrics{1});
        
        for i_metric_index = 1:length(metrics)
            for i_metric_name = 1:length(metric_names)
                metrics_cal{i_metric_index}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
            end
        end
        
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        % codon_list = {'TotalActivity'};
        
        
        if single_cell_specificity_plot
            clear codon_single_cell
            pairwise_ligand_index = {[1,2];[1,3];[1,4];[1,5];[2,3];[2,4];[2,5];[3,4];[3,5];[4,5]};
            codon_single_cell=cell(2,length(pairwise_ligand_index));
            
            for j_cell = 1:length(pairwise_ligand_index)
                for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                    for i_ligand = 1:length(pairwise_ligand_index{j_cell})
                        for i_codon = 1:length(codon_list)
                            codon_single_cell{i_cell,j_cell}(i_ligand,i_codon) = collect_feature_vects.(codon_list{i_codon}){pairwise_ligand_index{j_cell}(i_ligand)}(i_cell,:);
                        end
                    end
                end
            end
            
            clusterResults = cell2mat(cellfun(@applyDBSCAN, codon_single_cell, 'UniformOutput', false));
            
            MI_bars = mean(log2(clusterResults));
            bar(MI_bars)
        end
        
        
    end
end

%% epsilon network clustering
function cluster_num = applyDBSCAN(matrix)
epsilon = 1; % Set the epsilon value
minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end