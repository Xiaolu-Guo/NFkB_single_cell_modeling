function    CpG_polyIC_compete_compare_supriya_202309(data_file)
% For Guo et al. Figure 4G, signaling codon distribution for sampling data
% with CpG-polyIC competition
% 
% tested 05/12/2024, Matlab 2020a

% clear all
% load the codon for all sampled data and supriya and ade's data for
% comparison
% comapring ade's data and supriya's data: the batch effects.
% comparing surpiya's data and model prediction
% check each chunk to set 1 or 0 to run or not run the corresponding part
%
cal_codon_all = 0;
cal_codon_dual_ligand = 1;
cal_codon = 1;


addpath('./lib/')
addpath('./src/')
addpath('./bin/')

fig_save_path = '../SubFigures2023/';

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';


%%  calculate/load the codon for dual ligand
if cal_codon_dual_ligand
    if cal_codon
        
        
        % sampled data
        data_filename = 'Sim5_codon_all5dose_metric.mat';
        load(strcat(data_save_file_path_1,data_filename))
        sample_index = [11;15];
        
        metric_fields_name = fieldnames(metrics{1});
        for i_metric = 1:length(metrics)
            
            for i_metric_fields = 1:length(metric_fields_name)
                metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
            end
        end
        
        data_sample = data;
        
        load(strcat(data_save_file_path_1,data_file))
        sample_dual_ligand_index = 1;
        
        metric_fields_name = fieldnames(metrics{1});
        for i_metric = 1:length(metrics)
            
            for i_metric_fields = 1:length(metric_fields_name)
                metrics_sample_dual{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
            end
        end
        
        data_sample_dual = data;
        
        % metrics_sample = metrics;
        
        metrics = [metrics_sample(sample_index),metrics_sample_dual(sample_dual_ligand_index)]; % sampling dual ligand data
        
        
        data_info.info_ligand = [data_sample.info_ligand(sample_index),data_sample.info_ligand(sample_dual_ligand_index)]; % sampling dual ligand data
        
        
        for i_info_ligand = 1:length(data_info.info_ligand)
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'polyIC','PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');
            data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_LPS','LPS_PolyIC');
        end
        
        
        data_info.info_dose_str = [data_sample.info_dose_str(sample_index),data_sample_dual.info_dose_str(sample_dual_ligand_index)]; % sampling dual ligand data
        
        data_info.data_label(1:length(sample_index)) = {'Sampling'};
        data_info.data_label(end+1:end+length(sample_dual_ligand_index)) = {'Sampling_dual_ligand'};
        
        % [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023(data_info,metrics); %,  parameter
        
        [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics); %,  parameter
        metrics = metrics_new;
        %
        collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
        collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
        collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
        %
        % save(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'),'collect_feature_vects','metrics','data_info')
        % save(strcat(data_save_file_path_1,'All_codon_dual_ligand_202307_nonminmaxscaled.mat'),'collect_feature_vects','metrics','data_info')
        
        % goal: visualize supriya, ade, fitting, sampling (2 seprate plots) as
        % benchmarking; then visualize supriya dual ligand vs sampling.
        
    end
end


%% violin plot of codon (duration and total activity) distribution under co stim vs single stim for comparison, only picked three
% Figure 7 v08302023

% A + B : [supriya A, supriya B, supriya A + B],  [sampling A, sampling B, sampling A + B]
% TNF + Pam3CSK : [1,2,11],  [6,7,20]
% TNF + CpG : [1,3,12],  [6,8,21]
% TNF + LPS : [1,4,13],  [6,9,22]
% TNF + PolyIC : [1,5,14],  [6,10,23]
% Pam3CSK + CpG : x,  [7,8,24]
% Pam3CSK + LPS : [2,4,15],  [7,9,25]
% Pam3CSK + PolyIC : [2,5,16],  [7,10,26]
% CpG + LPS : [3,4,17],  [8,9,27]
% CpG + PolyIC : [3,5,18],  [8,10,28]
% LPS + PolyIC : [4,5,19],  [9,10,29]

if 1
    codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
    
    clear index_violin collect_feature_vects_draw
    codon_fields = fieldnames(collect_feature_vects);
    
    index_violin_sampling = {[1,2,3]};%pick
    
    
    index_violin_vec= {index_violin_sampling};
    color_vec = {[0,1,0]};
    data_type = {'sampling'};
    
    for i_index = 1:length(index_violin_vec)
        index_violin = index_violin_vec{i_index};
        for i_index_violin = 1:length(index_violin)
            clear collect_feature_vects_draw
            
            for i_codon_fields = 1:length(codon_fields)
                collect_feature_vects_draw.(codon_fields{i_codon_fields}) = collect_feature_vects.(codon_fields{i_codon_fields})(index_violin{i_index_violin});
            end
            
            for i_codon = 1:length(codon_list)
                collect_feature_vects_draw.(codon_list{i_codon}) = combined_cell_zscore(collect_feature_vects_draw.(codon_list{i_codon}));
            end
            
            
            for i_data_types = 1:length(collect_feature_vects_draw.info_data_type)
                collect_feature_vects_draw.info_data_type{i_data_types} = data_type{i_index};
            end
            
            fig_opt.save_file = strcat(fig_save_path,'CpG_polyIC_competetion_',data_type{i_index},'_compare_codon_distrib_202309');
            %     fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
            %     fig_opt.paper_opt.papersize=[220 180]*3;
            fig_opt.paper_opt.paperpos=[0,0,100,50];
            fig_opt.paper_opt.papersize=[100 50];
            fig_opt.distri_color = color_vec(i_index);
            fig_opt.codon = {'Duration','TotalActivity'};
            
            % fig_opt.save_file = strcat(fig_save_path,'Co_Stim_compare_codon_distrib_202306');
            if 1
                    violin_plot_codon_compare_co_sti_202308(collect_feature_vects_draw,fig_opt)
  
            end
            
            if 1
                codon_list = {'Duration','TotalActivity'};
                
                for i_codon=1:length(codon_list)
                    
                    vects = collect_feature_vects_draw.(codon_list{i_codon});
                    median_codon = cellfun(@median,vects(1:3));
                    codon_list{i_codon}
                    collect_feature_vects_draw.info_ligand{3}
                    collect_feature_vects_draw.info_data_type{3}
                    diff = median_codon(3) - max(median_codon(1:2))
                    
                end
            end
            
        end
    end
    
end


end


%% funciton zscore

function zscores = combined_cell_zscore(dataCells)
% Concatenate all data from the cells to compute the overall mean and standard deviation
allData = vertcat(dataCells{:});
mu = nanmean(allData);
sigma = nanstd(allData);

% Initialize the output cell array with the same size as the input
zscores = cell(size(dataCells));

% Loop through each cell and compute the z-scores for the data inside using the overall mean and std
for i = 1:numel(dataCells)
    zscores{i} = (dataCells{i} - mu) / sigma;
end
end
