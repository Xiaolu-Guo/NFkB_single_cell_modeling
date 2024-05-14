clear all

%% initializing do not copy this
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';
savepath = fig_save_path;
addpath('./lib/')
addpath('./src/')
addpath('./bin/')

%% draw 5 ligand single ligand stim

vers = '_r1';
vers_savefig = strcat(vers,'_matching_SS_p25x');%_IkBao _IkBao_onlyPeak ,'_SS_0x' _p1x
if 1 %draw 5 ligand single ligand stim
    % Sim16_IkBao_5_signle_ligand_codon_metric _r1?
    % Sim15_5_signle_ligand_codon_metric _r2
    % Sim16_IkBao_5_signle_ligand_codon_metric_0x r1
    % Sim16_IkBao_5_signle_ligand_codon_metric_p01x
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_0x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p1x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x _r1
    % Sim16_IkBao_matching_5_signle_ligand_codon_metric_p01x _r1
    load(strcat('../raw_data2023/Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x',vers,'.mat'))
    data_IkBo = data;
    para_mat_IkBo = sim_data_tbl.parameter_value(1:9:end,:);
    %% draw codons
    if 1
        % codon_mat = double;
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
        
        length_data = length(data.model_sim);
        length_metrics = length(metrics);
        
        load(strcat('../raw_data2023/Sim8_5_signle_ligand_codon_metric_r3','.mat'))
        para_mat_wt = sim_data_tbl.parameter_value(1:9:end,:);
        
        for i_data = 1:length(data.model_sim)
            
            data_info.info_ligand{i_data+length_data} = data.info_ligand{i_data};
            data_info.info_dose_str{i_data+length_data} = data.info_dose_str{i_data};
            data_info.data_label{i_data+length_data} = 'SamplingComb';
            
        end
        
        for i_metric_index = 1:length(metrics)
            for i_metric_name = 1:length(metric_names)
                metrics_cal{i_metric_index+length_metrics}.(metric_names{i_metric_name}) = metrics{i_metric_index}.(metric_names{i_metric_name})(1:9:end,:);
            end
        end
        
        
        %[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        [collect_feature_vects,metrics] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %calculate_codon_from_metric2023_07_nonminmaxscaled
        
        codon_list = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};
        weights_codon = [1,1,1,1,1,1];
        %codon_list = {'TotalActivity'};
        
        %% draw single cells
        if single_cell_specificity_plot
            clear codon_single_cell_wt codon_single_cell_SS
            codon_single_cell_wt=cell(2,1);
            codon_single_cell_IkBo=cell(2,1);
            
            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:length_data
                    for i_codon = 1:length(codon_list)
                        codon_single_cell_IkBo{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand}(i_cell,:);
                    end
                end
            end
            
            for i_cell = 1:size(collect_feature_vects.(codon_list{1}){1},1)
                for i_ligand = 1:length_data
                    for i_codon = 1:length(codon_list)
                        codon_single_cell_wt{i_cell}(i_ligand,i_codon) = weights_codon(i_codon)*collect_feature_vects.(codon_list{i_codon}){i_ligand+length_data}(i_cell,:);
                    end
                end
            end
            
            
            if 1
                dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
                dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
                
                % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
                %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
                %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
                %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
                
                %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
                %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
                %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
                
                dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
                dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
                
                
                
                if 0 
                figure(1)
                paperpos=[0,0,250,100]*1.5;
                papersize=[250 100]*1.5;
                draw_pos=[10,10,230,90]*1.5;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                
                y = dist_mat_wt;
                z = dist_mat_IkBo;
                % subplot(1,length(vis_data_field),i_data_field)
                
                al_goodplot_pair_RMSD(y,[],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'left',[],std(y(:))/2500);
                al_goodplot_pair_RMSD(z,[],0.5,ones(size(z,2),1)*[255 0 0]/255,'right',[],std(z(:))/2500);
                
                xlim([0.4 10.6])
                
                xticks([1:10])
                xticklabels({})
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
                
                ylim([0,5]);
                for i_x = 1:10
                    plot([i_x,i_x],[0,5],'--','Color','k');hold on
                end
                set(gca,'fontsize',14,'fontname','Arial');
                % saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');
                saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_10pert_med',vers_savefig),'epsc');

                close
                
                end
                
                if 0
                figure(1)
                paperpos=[0,0,250,100]*1.5;
                papersize=[250 100]*1.5;
                draw_pos=[10,10,230,90]*1.5;
                set(gcf, 'PaperUnits','points')
                set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                
                y = dist_mat_wt;
                z = dist_mat_wt;
                % subplot(1,length(vis_data_field),i_data_field)
                
                al_goodplot_pair_RMSD(y,[],0.5,ones(size(y,2),1)*[0 0 255]/255 ,'left',[],std(y(:))/2500);
                al_goodplot_pair_RMSD(z,[],0.5,ones(size(z,2),1)*[0 0 255]/255,'right',[],std(z(:))/2500);
                
                xlim([0.4 10.6])
                
                xticks([1:10])
                xticklabels({})
                %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
                
                ylim([0,5]);
                for i_x = 1:10
                    plot([i_x,i_x],[0,5],'--','Color','k');hold on
                end
                set(gca,'fontsize',14,'fontname','Arial');
                % saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');
                saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_wt_exp_10pert_med',vers_savefig),'epsc');

                close
                
                end
            end
            
            if 1
                dist_mat_wt = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_wt, 'UniformOutput', false));
                dist_mat_IkBo = cell2mat(cellfun(@(x) dist_mat(x), codon_single_cell_IkBo, 'UniformOutput', false));
                
                % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
                %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
                %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
                %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
                
                conf_pairs = {'TNF-LPS', 'TNF-CpG', 'TNF-P3C', 'TNF-PIC',... [1,2,4,3]
                    'LPS-CpG', 'LPS-P3C', 'CpG-P3C',... [5,7,9]
                    'LPS-PIC', 'CpG-PIC', 'P3C-PIC'};% [6,8,10]
                
                %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
                %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
                %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
                
                dist_mat_wt = dist_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
                dist_mat_IkBo = dist_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
                % 805, RMSD approx 0.5; 550, RMSD approx 1.5
                
                if 1
                    
                    [fig1,fig2,fig3,conf_pairs] = hierarchical_row_column_plot(dist_mat_wt,[0,5],parula,conf_pairs);
                    
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_wt_hierarchical_cluster_rows'),'epsc');%_20cells
                    close
                    
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_wt_hierarchical_cluster_columns'),'epsc');%_20cells
                    close
                    
                    fig3
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_wt_hierarchical_heatmaps'),'epsc');%_20cells
                    close
                    
                    [fig1,fig2,fig3] = hierarchical_row_column_plot(dist_mat_IkBo,[0,5],parula,conf_pairs,0);
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_IkBs_hierarchical_cluster_rows'),'epsc');%_20cells
                    close
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_IkBs_hierarchical_cluster_columns'),'epsc');%_20cells
                    close
                    fig3
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_IkBs_hierarchical_heatmaps'),'epsc');%_20cells
                    close
                    
                end
                
                if 0 % representative for low and high l2 distance
                    figure(1)
                    paperpos=[0,0,120,100];
                    papersize=[120 100];
                    draw_pos=[10,10,100,90];
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    plot(0:1/12:8,metrics{6}.time_series(805,:),'b','LineWidth',1.5); hold on
                    plot(0:1/12:8,metrics{7}.time_series(805,:),'r','LineWidth',1.5)
                    xticks(0:2:8)
                    xlim([0,8])
                    xticklabels('')
                    yticks(0:0.1:0.4)
                    ylim([0,0.4])
                    yticklabels('')
                    saveas(gcf,strcat(savepath,'representative_low_Euclidean'),'epsc');%_20cells
                    close
                    
                    figure(2)
                    paperpos=[0,0,120,100];
                    papersize=[120 100];
                    draw_pos=[10,10,100,90];
                    set(gcf, 'PaperUnits','points')
                    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                    plot(0:1/12:8,metrics{6}.time_series(578,:),'b','LineWidth',1.5); hold on
                    plot(0:1/12:8,metrics{7}.time_series(578,:),'r','LineWidth',1.5)
                    xticks(0:2:8)
                    xlim([0,8])
                    ylim([0,0.4])
                    xticklabels('')
                    yticks(0:0.1:0.4)
                    yticklabels('')
                    saveas(gcf,strcat(savepath,'representative_high_Euclidean'),'epsc');%_20cells
                    close
                end
                
                if 0
                    
                    [fig1,fig2] = hierarchical_plots(dist_mat_wt);
                    
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_wt_hierarchical_cluster'),'epsc');%_20cells
                    close
                    
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_wt_heatmaps'),'epsc');%_20cells
                    close
                    
                    
                    [fig1,fig2] = hierarchical_plots(dist_mat_IkBo);
                    
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_IkBs_hierarchical_cluster'),'epsc');%_20cells
                    close
                    
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_RMSD_IkBs_heatmaps'),'epsc');%_20cells
                    close
                end
                
                
                if 0
                    RMSD_threshold = 1;
                    conf_mat_wt = cell2mat(cellfun(@(x) confusion_mat(x,RMSD_threshold), codon_single_cell_wt, 'UniformOutput', false));
                    [fig1,fig2] = hierarchical_plots(conf_mat_wt,[0,1],[1,1,1;0,0,0]);
                    
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_wt_hierarchical_cluster'),'epsc');%_20cells
                    close
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_wt_heatmaps'),'epsc');%_20cells
                    close
                    
                    
                    conf_mat_IkBo = cell2mat(cellfun(@(x) confusion_mat(x,RMSD_threshold), codon_single_cell_IkBo, 'UniformOutput', false));
                    [fig1,fig2] = hierarchical_plots(conf_mat_IkBo,[0,1],[1,1,1;0,0,0]);
                    
                    fig2
                    saveas(gcf,strcat(savepath,'Confusion_pairs_IkBs_hierarchical_cluster'),'epsc');%_20cells
                    close
                    fig1
                    saveas(gcf,strcat(savepath,'Confusion_pairs_IkBs_heatmaps'),'epsc');%_20cells
                    close
                end
                
                
            end
            
            %%
            if 0
                for epsilon = 1 % :3
                    % epsilon = 2; % Define the epsilon value
                    
                    conf_mat_wt = cell2mat(cellfun(@(x) confusion_mat(x,epsilon), codon_single_cell_wt, 'UniformOutput', false));
                    conf_mat_IkBo = cell2mat(cellfun(@(x) confusion_mat(x,epsilon), codon_single_cell_IkBo, 'UniformOutput', false));
                    
                    % {'TNF';'LPS';'CpG';'PolyIC';'Pam3CSK'}
                    %TNF-LPS, TNF-CpG, TNF-PIC ,TNF-P3C [1,2,3,4]
                    %LPS-CpG, LPS-PIC, LPS-P3C, [5,6,7]
                    %CpG-PIC, CpG-P3C, P3C-PIC,[8,9,10]
                    
                    %TNF-LPS, TNF-CpG, TNF-P3C, TNF-PIC [1,2,4,3]
                    %LPS-CpG, LPS-P3C, CpG-P3C [5,7,9]
                    %LPS-PIC, CpG-PIC, P3C-PIC [6,8,10]
                    
                    conf_mat_wt = conf_mat_wt(:,[1,2,4,3,5,7,9,6,8,10]);
                    conf_mat_IkBo = conf_mat_IkBo(:,[1,2,4,3,5,7,9,6,8,10]);
                    [ax,fig] = plot_bar_conf(conf_mat_wt,size(conf_mat_wt,1),[]);
                    if epsilon ==1
                        ylim([0 0.3]);
                        set(gca,'YTick',0:0.15:0.3)%'XTick',0:1:5,
                    else
                        ylim([0 1]);
                        set(gca,'YTick',0:0.5:1)%'XTick',0:1:5,
                    end
                    saveas(gcf,strcat(fig_save_path,'SRS_conf_counts',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                    close
                    
                    if 1
                        conf_ligand_wt = sum(conf_mat_wt)/size(conf_mat_wt,1);
                        conf_ligand_IkBo = sum(conf_mat_IkBo)/size(conf_mat_IkBo,1);
                        
                        figure(1)
                        paperpos = [0,0,100,50]*2;
                        papersize = [100,50]*2;
                        draw_pos = [10,10,90,30]*2;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        bar([conf_ligand_wt;conf_ligand_IkBo]')
                        
                        ylim([0 0.25]);
                        set(gca,'YTick',0:0.05:0.25)%'XTick',0:1:5,
                        
                        % Set x-axis tick labels to none
                        xticklabels({});
                        % Set y-axis tick labels to none
                        yticklabels({});
                        
                        if epsilon ==1
                            ylim([0 0.3]);
                            set(gca,'YTick',0:0.15:0.3)%'XTick',0:1:5,
                        else
                            ylim([0 1]);
                            set(gca,'YTick',0:0.5:1)%'XTick',0:1:5,
                        end
                        saveas(gcf,strcat(fig_save_path,'SRS_conf_counts_IkBo_wt_pairs_',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                    
                    conf_score_wt = sum(conf_mat_wt,2);
                    conf_score_IkBo = sum(conf_mat_IkBo,2);
                    
                    for cluster_num =0:10
                        barplot(cluster_num+1,1) = sum(conf_score_wt == cluster_num) ;
                        barplot(cluster_num+1,2) = sum(conf_score_IkBo == cluster_num) ;
                    end
                    
                    % epsilon
                    % barplot
                    
                    if 0
                        
                        figure(1)
                        paperpos = [0,0,80,50]*1.5;
                        papersize = [80,50]*1.5;
                        draw_pos = [10,10,60,30]*1.5;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        bar(barplot(:,1))
                        
                        ylim([0,500])
                        set(gca,'YTick',[0,250,500])
                        
                        % Set x-axis tick labels to none
                        xticklabels({});
                        
                        % Set y-axis tick labels to none
                        yticklabels({});
                        
                        saveas(gcf,strcat(fig_save_path,'Single_cell_conf_score_bar_wt',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                    
                    if 1
                        
                        figure(1)
                        paperpos = [0,0,80,50]*1.5;
                        papersize = [80,50]*1.5;
                        draw_pos = [10,10,60,30]*1.5;
                        set(gcf, 'PaperUnits','points')
                        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
                        
                        bar(barplot)
                        
                        
                        
                        ylim([0,500])
                        set(gca,'YTick',[0,250,500])
                        % Set x-axis tick labels to none
                        xticklabels({});
                        
                        % Set y-axis tick labels to none
                        yticklabels({});
                        
                        saveas(gcf,strcat(fig_save_path,'Single_cell_conf_score_bar_wt_IkBs',vers_savefig,'_eps',num2str(epsilon)),'epsc')
                        close
                        
                    end
                    
                end
                
            end
            
            
            
            
        end
    end
    
    %% draw heat map IkBo-/-
    if 0
        clear data_NFkB
        cells_inteval = 50;
        
               ligand_index = [1,5,3,2,4]; % reoder the stim

        data_NFkB.info_ligand = data_IkBo.info_ligand(ligand_index);
        data_NFkB.info_dose_str = data_IkBo.info_dose_str(ligand_index);
        data_NFkB.info_dose_index = data_IkBo.info_dose_index(ligand_index);
        
        for i_ligand = 1:length(data_IkBo.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data_IkBo.model_sim{ligand_index(i_ligand)}(1:9:end,:);
            
        end
        
        data_plot = [];
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_plot = [data_plot,data_NFkB.model_sim{i_ligand}(1:cells_inteval:end,:)];
            %data_NFkB.exp = data_NFkB.model_sim;
        end
        
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        [~,order_to_plot] = sort(metrics{1}.(metric_order_name_{1})(1:cells_inteval:end,metric_order_column),'descend');
        
        figure(2)
        
        clear traj_mat;
        traj_mat = data_plot;
        
        % Step 1: Perform hierarchical clustering on the rows
        Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
        Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
        % Step 2: Determine the order of rows based on hierarchical clustering
        [H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
        close(gcf); % Close dendrogram figure
        
        % Step 3: Reorder the matrix based on the clustering result
        reordered_traj_mat = traj_mat(Outperm, :);
        
        % If you need to display the dendrogram alongside, you can plot it separately
        figure(2)
        paperpos=[0,0,55,70];
        papersize=[55 70];
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        [dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
        set(gca, 'XDir', 'reverse','YDir','reverse');
        set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
        xticklabels('')
        yticklabels('')
        axis off
        saveas(gcf,strcat(savepath,'SRS_heatmap_IkBm_hierarchical_cluster_20cells'),'epsc');%_20cells
        close
        
        
        figure(1)
        
        paperpos=[0,0,246,64]*3;
        papersize=[246,64]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        %
        XLabels = 0:5:((size(reordered_traj_mat,2)-1)*5);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/60);
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:size(data_plot,1);
        % Convert each number in the array into a string
        YCustomXLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        YCustomXLabels(:) = " ";
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.YDisplayLabels = YCustomXLabels;
        
        % xlabel('Time (hours)');
        % ylabel(vis_data_field{i_data_field});
        % clb=colorbar;
        % clb.Label.String = 'NFkB(A.U.)';
        colorbar('off')
        
        set(gca,'fontsize',14,'fontname','Arial');
        saveas(gcf,strcat(savepath,'SRS_heatmap_IkBm_20cells'),'epsc');%_20cells
        close
        
    end
    
    if 0
        data_NFkB.info_ligand = data_IkBo.info_ligand;
        data_NFkB.info_dose_str = data_IkBo.info_dose_str;
        data_NFkB.info_dose_index = data_IkBo.info_dose_index;
        for i_ligand = 1:length(data_IkBo.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data_IkBo.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        
        order_name_vec = {'matched_sc_srs_TNF','matched_sc_srs_LPS','matched_sc_srs_CpG','matched_sc_srs_PolyIC','matched_sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers_savefig);
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 1 % 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{i_order_i_ligand}.(metric_order_name_{1})(:,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
    %% draw heat map wt
    if 0
        clear data_NFkB
        cells_inteval = 50;
       ligand_index = [1,5,3,2,4]; % reoder the stim

        data_NFkB.info_ligand = data.info_ligand(ligand_index);
        data_NFkB.info_dose_str = data.info_dose_str(ligand_index);
        data_NFkB.info_dose_index = data.info_dose_index(ligand_index);

        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data.model_sim{ligand_index(i_ligand)}(1:9:end,:);
            
        end
        data_plot = [];
        
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_plot = [data_plot,data_NFkB.model_sim{i_ligand}(1:cells_inteval:end,:)];
            %data_NFkB.exp = data_NFkB.model_sim;
        end
        
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        [~,order_to_plot] = sort(metrics{length_metrics+1}.(metric_order_name_{1})(1:cells_inteval:end,metric_order_column),'descend');
        
        figure(2)
        clear traj_mat;
        traj_mat = data_plot;
        
        % Step 1: Perform hierarchical clustering on the rows
        Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
        Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
        % Step 2: Determine the order of rows based on hierarchical clustering
        [H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
        close(gcf); % Close dendrogram figure
        
        % Step 3: Reorder the matrix based on the clustering result
        reordered_traj_mat = traj_mat(Outperm, :);
        
        % If you need to display the dendrogram alongside, you can plot it separately
        figure(2)
        paperpos=[0,0,55,70];
        papersize=[55 70];
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
        
        [dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
        set(gca, 'XDir', 'reverse','YDir','reverse');
        set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
        
        xticklabels('')
        yticklabels('')
        axis off
        saveas(gcf,strcat(savepath,'SRS_heatmap_wt_hierarchical_cluster_20cells'),'epsc');%_20cells
        close
        
        
        figure(1)
        
        paperpos=[0,0,246,64]*3;
        papersize=[246,64]*3;
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)
        
        % subplot(1,length(vis_data_field),i_data_field)
        h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[-0.001,0.25]);%[-0.001,0.2] for TNF
        %
        XLabels = 0:5:((size(reordered_traj_mat,2)-1)*5);
        % Convert each number in the array into a string
        CustomXLabels = string(XLabels/60);
        % Replace all but the fifth elements by spaces
        % CustomXLabels(mod(XLabels,60) ~= 0) = " ";
        CustomXLabels(:) = " ";
        
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.XDisplayLabels = CustomXLabels;
        
        YLabels = 1:size(reordered_traj_mat,1);
        % Convert each number in the array into a string
        YCustomXLabels = string(YLabels);
        % Replace all but the fifth elements by spaces
        YCustomXLabels(:) = " ";
        % Set the 'XDisplayLabels' property of the heatmap
        % object 'h' to the custom x-axis tick labels
        h.YDisplayLabels = YCustomXLabels;
        
        % xlabel('Time (hours)');
        % ylabel(vis_data_field{i_data_field});
        % clb=colorbar;
        % clb.Label.String = 'NFkB(A.U.)';
        colorbar('off')
        
        set(gca,'fontsize',14,'fontname','Arial');
        saveas(gcf,strcat(savepath,'SRS_heatmap_wt_20cells'),'epsc');%_20cells
        close
        
    end
    
    if 0
        clear data_NFkB
        data_NFkB.info_ligand = data.info_ligand;
        data_NFkB.info_dose_str = data.info_dose_str;
        data_NFkB.info_dose_index = data.info_dose_index;
        for i_ligand = 1:length(data_NFkB.info_ligand)
            data_NFkB.info_num_cells{i_ligand} = 1000;
            data_NFkB.model_sim{i_ligand} = data.model_sim{i_ligand}(1:9:end,:);
            
        end
        
        data_NFkB.exp = data_NFkB.model_sim;
        vis_data_field = {'model_sim'};
        % integrals 97 oscfreq 1 oscpower 1
        metric_order_name_{1} = 'integrals';
        metric_order_column = 97;
        
        vers_savefig = '_matching_wt';
        order_name_vec = {'matched_sc_srs_TNF','matched_sc_srs_LPS','matched_sc_srs_CpG','matched_sc_srs_PolyIC','matched_sc_srs_Pam3CSK'};
        for i_order_name_vec = 1:length(order_name_vec)
            codon_name_vec{i_order_name_vec} = strcat('_',metric_order_name_{1},'_',vers_savefig);
            
        end
        order_name_vec = cellfun(@strcat, order_name_vec,codon_name_vec,'UniformOutput',false);
        
        % 'sc_srs_' represents single-cell stumilus response specificity
        
        for i_order_i_ligand = 1% 1:length(data_NFkB.info_ligand)
            
            for i_ligand = 1:length(data_NFkB.info_ligand)
                [~,data_NFkB.order{i_ligand}] = sort(metrics{length_metrics+i_order_i_ligand}.(metric_order_name_{1})(:,metric_order_column),'descend');
            end
            filter_TNF = 0;
            plot_traj_heatmap_2024_samplingSRS(data_NFkB,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_i_ligand})
        end
        
    end
    
end

%%
function cluster_num = applyDBSCAN(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
cluster_num = max(idx);
end

%%
function conf_mat = confusion_mat(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3
conf_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        conf_mat(i_col) = (sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2))<epsilon);
        i_col = i_col + 1;
    end
end

end

%%
function dist_mat = dist_mat(matrix)

% epsilon = 2; % Set the epsilon value 1,2,3
dist_mat = -ones(1,size(matrix,1) * (size(matrix,1)-1)/2);
i_col = 1;
for i_sti1 = 1:size(matrix,1)
    for i_sti2 = (i_sti1+1):size(matrix,1)
        dist_mat(i_col) = sqrt(sum((matrix(i_sti1,:) - matrix(i_sti2,:)).^2));
        i_col = i_col + 1;
    end
end

end

%%
function cluster_results = applyDBSCAN_clusters(matrix,epsilon)

% epsilon = 2; % Set the epsilon value 1,2,3

minPts = 1;    % Set the minimum number of points
idx = dbscan(matrix, epsilon, minPts);
i_cluster = 1;
idx_rest = idx';
for i_cluster = 1:max(idx)
    cluster_results{i_cluster} = find(idx_rest == i_cluster);
end
end

%%
function confusion_mat = confusion_matrix_cal(clusters_genetype,all_combine)
confusion_mat = zeros(length(clusters_genetype),length(all_combine));
for i_cell = 1:length(clusters_genetype)
    for i_2combine = 1:length(all_combine)
        for i_clusters = 1:length(clusters_genetype{i_cell})
            if isequal(sort(all_combine{i_2combine}), sort(clusters_genetype{i_cell}{i_clusters}))
                confusion_mat(i_cell,i_2combine) = 1;
            end
            
        end
        
    end
    
end
end

%%
function confusion_mat = confusion_all_matrix_cal(clusters_genetype,all_combine)
confusion_mat = zeros(length(clusters_genetype),length(all_combine));
for i_cell = 1:length(clusters_genetype)
    for i_2combine = 1:length(all_combine)
        for i_clusters = 1:length(clusters_genetype{i_cell})
            if all(ismember(all_combine{i_2combine}, clusters_genetype{i_cell}{i_clusters}))
                confusion_mat(i_cell,i_2combine) = 1;
            end
        end
    end
end

end

%%
function [ax,fig] = plot_bar_conf(confusion_mat,total_cell,idx_cell)

if isempty(total_cell)
    total_cell = size(confusion_mat,1);
end

if isempty(idx_cell)
    conf_ligand = sum(confusion_mat)/total_cell;
    
else
    conf_ligand = sum(confusion_mat(idx_cell,:))/total_cell;
end

figure(1)
paperpos = [0,0,100,50]*2;
papersize = [100,50]*2;
draw_pos = [10,10,90,30]*2;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
conf_ligand
bar(conf_ligand)

ylim([0 0.25]);
set(gca,'YTick',0:0.05:0.25)%'XTick',0:1:5,

% Set x-axis tick labels to none
xticklabels({});

% Set y-axis tick labels to none
yticklabels({});
ax = gca;
fig = gcf;
end

%%
function [fig1,fig2] = hierarchical_plots(traj_mat,color_limit,color_map)

figure(1)

% Step 1: Perform hierarchical clustering on the rows
Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% Step 3: Reorder the matrix based on the clustering result
reordered_traj_mat = traj_mat(Outperm, :);

% If you need to display the dendrogram alongside, you can plot it separately
figure(1)
paperpos=[0,0,55,70];
papersize=[55 70];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig1 = gcf;

figure(2)

paperpos=[0,0,100,64]*3;
papersize=[100,64]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

% subplot(1,length(vis_data_field),i_data_field)
if nargin <2
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,5]);%[-0.001,0.2] for TNF
elseif nargin <3
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
else
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',color_map,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
    
end
%
XLabels = 1:size(reordered_traj_mat,2);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/1);
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_mat,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig2 = gcf;
end
%%
function [fig1,fig2,fig3,conf_pairs] = hierarchical_row_column_plot(traj_mat,color_limit,color_map,conf_pairs,order_conf_pairs)

figure(1)

% Step 1: Perform hierarchical clustering on the rows
Y = pdist(traj_mat, 'euclidean'); % Compute the pairwise distances between rows
Z = linkage(Y, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm] = dendrogram(Z, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(1)
paperpos=[0,0,55,70];
papersize=[55 70];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h,~,~] = dendrogram(Z,0, 'Orientation', 'left'); % Plotting dendrogram separately
set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig1 = gcf;

figure(2)
% Step 1: Perform hierarchical clustering on the rows
Y_col = pdist(traj_mat', 'euclidean'); % Compute the pairwise distances between rows
Z_col = linkage(Y_col, 'ward'); % Perform hierarchical/agglomerative clustering
% Step 2: Determine the order of rows based on hierarchical clustering
[H,T,Outperm_col] = dendrogram(Z_col, 0); % Get the order of rows for clustering
close(gcf); % Close dendrogram figure

% If you need to display the dendrogram alongside, you can plot it separately
figure(2)
paperpos=[0,0,70,55];
papersize=[70 55];
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

[dendro_h_col,~,~] = dendrogram(Z_col,0, 'Orientation', 'top'); % Plotting dendrogram separately
% set(gca, 'XDir', 'reverse','YDir','reverse');
set(dendro_h_col,'LineWidth',0.75,'Color','k'); % Adjust line width for better visibility
xticklabels('')
yticklabels('')
axis off
fig2 = gcf;

% Step 3: Reorder the matrix based on the clustering result
reordered_traj_mat = traj_mat(Outperm, Outperm_col);



figure(3)

paperpos=[0,0,100,64]*3;
papersize=[100,64]*3;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',paperpos)

% subplot(1,length(vis_data_field),i_data_field)
if nargin <2
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',[0,5]);%[-0.001,0.2] for TNF
elseif nargin <3
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',parula,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
else
    h=heatmap(reordered_traj_mat(:,:),'ColorMap',color_map,'GridVisible','off','ColorLimits',color_limit);%[-0.001,0.2] for TNF
    
end
%
XLabels = 1:size(reordered_traj_mat,2);
% Convert each number in the array into a string
CustomXLabels = string(XLabels/1);%conf_pairs
% Replace all but the fifth elements by spaces
% CustomXLabels(mod(XLabels,60) ~= 0) = " ";
CustomXLabels(:) = " ";
if nargin<5
    conf_pairs = conf_pairs(Outperm_col);
elseif order_conf_pairs
    conf_pairs = conf_pairs(Outperm_col);
end

CustomXLabels = string(conf_pairs);

% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_mat,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h.YDisplayLabels = YCustomXLabels;

% xlabel('Time (hours)');
% ylabel(vis_data_field{i_data_field});
% clb=colorbar;
% clb.Label.String = 'NFkB(A.U.)';
colorbar('off')

set(gca,'fontsize',14,'fontname','Arial');
fig3 = gcf;


end