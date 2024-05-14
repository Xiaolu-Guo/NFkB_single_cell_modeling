function [] =  save_sampling_20doses_for_channel_capacity(data_filename)

%This is for visualizing Monolix estimation

%% filepath


% clear all


% data_save_file_path = '../SAEM_proj_2022_2/';
% data_proj_num_vec = 5;

% load('allsti_data_codon_0603.mat')

%%

sample_data = 1; % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
vis_data_field = {'pred_mode_amp'};% 'pred_mode_filter_nan'};% 'model_sim_2'};%,'sample'};
vis_cv_field = {'vis_cv'};

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

% load(strcat(data_save_file_path_1,'Sim2_codon_metric.mat'))
%
% load(strcat(data_save_file_path_1,'Sim2_polyIC_codon_metric.mat'))
%
% load(strcat(data_save_file_path_1,'Sim2_r2_codon_metric.mat'))

% load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));

% load(strcat(data_save_file_path_1,'Sim2_fitting_lowdose_codon_metric.mat'))

% load(strcat(data_save_file_path_1,'Sim2_fitting_alldose_codon_metric.mat'))

load(strcat(data_save_file_path_1,data_filename))


% [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
%
% if sample_data % for sampled data, remove all the non NFkB trajectories, such as IKK traj.
%     data_fields_names = {'model_sim'};
%
%     for i_data = 1:length(data.(data_fields_names{1}))
%         for i_data_name = 1:length(data_fields_names)
%             data.(vis_data_field{1}){i_data} = data.(data_fields_names{i_data_name}){i_data}(1:9:end,:);
%         end
%     end
% else
%     data_fields_names = {'pred_mode_filter_nan'};
%     i_data_name =1;
%     data.(vis_data_field{i_data_name}) = data.(data_fields_names{i_data_name});
%
% end

metrics_new = cell(1,2);
metric_names = fieldnames(metrics{1});
for i_metric_name = 1:length(metric_names)
    for i_metric = 1:length(metrics)
        
        metrics_new{i_metric}.(metric_names{i_metric_name}) = metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
        data_info.info_ligand{i_metric} = collect_feature_vects.info_ligand{i_metric};
        data_info.info_dose_str{i_metric} = collect_feature_vects.info_dose_str{i_metric};
        data_info.data_label{i_metric} = 'sampling';
    end
end

[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_new); %,  parameter
collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);


%% save the signaling coodn data to mutual info format
clear nfkb
load('mutual_info_cal_data_example.mat')

codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};


sti_all = cellfun(@strcat,collect_feature_vects.info_ligand,collect_feature_vects.info_dose_str,'UniformOutput', false) ;
ligand_inf = collect_feature_vects.info_ligand{1};

for i_sti = 1:length(sti_all)
    nfkb(i_sti).sc_metrics = struct();
    
    for i_codon =1:length(codon_list)
        nfkb(i_sti).id = sti_all{i_sti};
        nfkb(i_sti).ids = sti_all;
        nfkb(i_sti).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){i_sti};
        
    end
end

% save(strcat('mutual_info_format_codon_',ligand_inf,'_20doses_20230614.mat'),'nfkb')

if 0
    switch ligand_inf
        case 'TNF'
            index_3_doses =[7,11,15];% d7,d11,d15
        case 'CpG' %: d9,d11,d13
            index_3_doses = [9,11,13];
        case 'LPS' %: d11,d13,d15
            index_3_doses = [11,13,15];
        case 'Pam3CSK' %: d7,d11,d15
            index_3_doses = [7,11,15];
        case 'PolyIC' %: d11,d13,d15
            index_3_doses = [11,13,15];
    end
    
    for i_index_3_doses = 1:length(index_3_doses)
        nfkb_3doses(i_index_3_doses) = nfkb(index_3_doses(i_index_3_doses));
        nfkb_3doses(i_index_3_doses).ids = nfkb(index_3_doses(i_index_3_doses)).ids(index_3_doses);
    end
    
    nfkb = nfkb_3doses;
    save(strcat('mutual_info_format_codon_',ligand_inf,'_3doses_20230614.mat'),'nfkb')
    
end

% save(strcat('mutual_info_format_codon_',ligand_inf,'_20doses_20230614.mat'),'nfkb')

if 1
    switch ligand_inf
        case 'TNF'
            index_3_doses =[1,7,11,15];% d7,d11,d15
        case 'CpG' %: d9,d11,d13
            index_3_doses = [1,9,11,13];
        case 'LPS' %: d11,d13,d15
            index_3_doses = [1,11,13,15];
        case 'Pam3CSK' %: d7,d11,d15
            index_3_doses = [1,7,11,15];
        case 'PolyIC' %: d11,d13,d15
            index_3_doses = [1,11,13,15];
    end
    
    for i_index_3_doses = 1:length(index_3_doses)
        nfkb_4doses(i_index_3_doses) = nfkb(index_3_doses(i_index_3_doses));
        nfkb_4doses(i_index_3_doses).ids = nfkb(index_3_doses(i_index_3_doses)).ids(index_3_doses);
    end
    
    nfkb = nfkb_4doses;
    save(strcat('mutual_info_format_codon_',ligand_inf,'_3doses_zerodose_20230707.mat'),'nfkb')
    
end

if 0
    index_3_doses = [1,11,21];
    for i_index_3_doses = 1:length(index_3_doses)
        nfkb_3doses(i_index_3_doses) = nfkb(index_3_doses(i_index_3_doses));
        nfkb_3doses(i_index_3_doses).ids = nfkb(index_3_doses(i_index_3_doses)).ids(index_3_doses);
    end
    
    nfkb = nfkb_3doses;
    save(strcat('mutual_info_format_codon_',ligand_inf,'_distinct_3doses_20230614.mat'),'nfkb')
    
end

