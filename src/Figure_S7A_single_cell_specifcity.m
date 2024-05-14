clear all

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

% Linearly interpolate between blue (at the start), white (in the middle), and red (at the end)
n =20;
blueToWhite = [linspace(0, 1, n/2)' linspace(0, 1, n/2)' ones(n/2, 1)];
whiteToRed = [ones(n/2, 1) linspace(1, 0, n/2)' linspace(1, 0, n/2)'];
customColormap = [blueToWhite; whiteToRed(2:end,:)];


%% draw 5 ligand single ligand stim

if 1 %draw 5 ligand single ligand stim
    load('../raw_data2023/Sim8_5_signle_ligand_codon_metric.mat')
    
    %% draw codons
    if 1
        % codon_mat = double;
        kmeans_cluster_plots = 1;
        single_cell_specificity_plot = 0;
        
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
        
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        % codon_list = {'TotalActivity'};
        
        
        %% draw kmean cluster
        if kmeans_cluster_plots
            
            clear codon_mat
            i_column = 1;
            ligand_index = [1,5,3,2,4]; % reoder the stim

            for i_ligand = 1:length(collect_feature_vects.info_ligand)
                for i_codon = 1:length(codon_list)
                    codon_mat(:,i_column) = collect_feature_vects.(codon_list{i_codon}){ligand_index(i_ligand)}(:,:);
                    i_column = i_column + 1;
                end
            end
            
            if 1
                % Assuming A is your 1000x30 matrix
                [rho,p] = corr(codon_mat, 'Type', 'Spearman');
                figure(1)
                h = heatmap(rho,'Colormap',customColormap);
                caxis([-1,1])
                
                YLabels = 1:size(codon_mat,2);
                % Convert each number in the array into a string
                YCustomXLabels = string(YLabels);
                % Replace all but the fifth elements by spaces
                YCustomXLabels(:) = " ";

                % Set the 'XDisplayLabels' property of the heatmap
                % object 'h' to the custom x-axis tick labels
                h.YDisplayLabels = YCustomXLabels;
                h.XDisplayLabels = YCustomXLabels;

                saveas(gcf,strcat(fig_save_path,'SRS_codon_corr_wt'),'epsc')
                close()
%                 figure(2)
%                 a = p<0.05;
%                 h = heatmap(double(a),'Colormap',customColormap);
%                 caxis([-1,1])
%                 
%                 YLabels = 1:size(codon_mat,2);
%                 % Convert each number in the array into a string
%                 YCustomXLabels = string(YLabels);
%                 % Replace all but the fifth elements by spaces
%                 YCustomXLabels(:) = " ";
% 
%                 % Set the 'XDisplayLabels' property of the heatmap
%                 % object 'h' to the custom x-axis tick labels
%                 h.YDisplayLabels = YCustomXLabels;
%                 h.XDisplayLabels = YCustomXLabels;
                
                a = 1;
             end
            
            if 0
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
        end
        %% draw single cells
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
            
            %             cell_idx = 58;
            %             figure(1)
            %             for i_ligand = 1:5
            %                 plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
            %             end
            %             ylabel('nNFkB (\mu M)')
            %             xlabel('Time (5min)')
            %             ylim([0,0.4])
            %             legend(collect_feature_vects.info_ligand)
            %             set(gca,'fontweight','b')
            %             %close()
            %             %specific traj
            %             % espilon = 2 => 709
            %
            %             % espilon = 3 => 709
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
            
            idx_vec = [668,58,709,752];
            mutual_info_vec = [0,1.59,2.32,2.32];
            
            if 0
                for i_idx = 1:4
                    cell_idx = idx_vec(i_idx);
                    cell_mutual_info = mutual_info_vec(i_idx);
                    figure(i_idx)
                    paperpos = [0,0,70,50]*1.8;
                    papersize = [70,50]*1.8;
                    draw_pos = [10,10,50,30]*1.8;
                    
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    
                    for i_ligand = 1:5
                        plot(1:length(metrics_cal{i_ligand}.time_series(cell_idx,:)),metrics_cal{i_ligand}.time_series(cell_idx,:),'LineWidth',2);hold on
                    end
                    title(strcat('cell MI = ',num2str(cell_mutual_info),'bits'))
                    ylabel('nNFkB (\mu M)')
                    xlabel('Time (5min)')
                    %lgd = legend(collect_feature_vects.info_ligand);
                    set(gca,'fontsize',6,'fontweight','b')
                    % Adding a legend
                    % Setting the legend box to be off
                    %set(lgd, 'Box', 'off');
                    saveas(gcf,strcat(fig_save_path,'Single_cell_spec_idx_',num2str(cell_idx)),'epsc')
                    close
                end
            end
            
            figure(1)
            histogram(clusterResults,5)
            
            % Calculating the mean
            meanValue = mean(clusterResults);
            
            
            % Adding a line at the mean value
            hold on; % Retain the current plot when adding new plots
            xline(meanValue, 'r', 'LineWidth', 2); % Red line for the mean value
            xlabel('cluster nums')
            ylabel('counts')
            
            figure(2)
            histogram(log2(clusterResults),5)
            meanValue = mean(log2(clusterResults))
            hold on; % Retain the current plot when adding new plots
            xline(meanValue, 'r', 'LineWidth', 2); % Red line for the mean value
            xlabel('cell MI (bits)')
            ylabel('counts')
            
            
            
            index = [9;9;9;9;9]*idx_vec + diag([0,9000,18000,27000,36000])*ones(5,4);
            index_vec = index(:);
            params_val = sim_data_tbl.parameter_value(index_vec,:);
            params_react_num = sim_data_tbl.parameter_reac_num(index_vec,:);
            params_para_num = sim_data_tbl.parameter_para_num(index_vec,:);
            
            params_val = params_val(1:5:end,:);
            params_react_num = params_react_num(1:5:end,:);
            params_para_num = params_para_num(1:5:end,:);
            
            
            
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

%%
function cluster_num = applyDBSCAN(matrix)
epsilon = 0.7; % Set the epsilon value
minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end