

data_exp =readmatrix(strcat('/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/Experiments/202202_rescaled_byXiaolu/','non_stim.xls'));
clear data
i_sti = 1;
data.exp{i_sti} = data_exp;
data.info_ligand{i_sti} = 'none';
data.info_dose_str{i_sti} = '0';
% cd(codon_info.codon_path)
vis_data_field = {'exp'};%,'sample'};
data_label = {'experiments'};%,'sample'};
[collect_feature_vects_unstim,metrics_unstim] = calculate_codon(data,vis_data_field,data_label);%,  parameter

fitting_dose ='alldose' ;%'lowdose','highdose','alldose','middose'
fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
fig_opt.paper_opt.papersize=[220 180]*3;

% data_proj_num_vec = data_proj_nums{i_ligand};
%     load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
% load(strcat(data_save_file_path,'Sim2_codon_metric.mat'))
data_save_file_path = '../raw_data2023/';%_fay_parameter/';
load(strcat(data_save_file_path,'Sim_unstim_fitting_alldose_r1_codon_metric.mat'))

metric_names = fieldnames(metrics{1});

for i_metric_name = 1:length(metric_names)
    metics_sim_unstim{1}.(metric_names{i_metric_name}) = [];
    for i_metric = 1:length(metrics)
        
        metics_sim_unstim{1}.(metric_names{i_metric_name}) = [metics_sim_unstim{1}.(metric_names{i_metric_name});metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:)];
    end
end


load(strcat(data_save_file_path,'Sim2_fitting_',fitting_dose,'_codon_metric.mat'))

% 1. read data and metric
% asign data
metric_names = fieldnames(metrics{1});

for i_metric = 1:length(metrics)
    for i_metric_name = 1:length(metric_names)
        metrics_new{i_metric}.(metric_names{i_metric_name}) = metrics{i_metric}.(metric_names{i_metric_name})(1:9:end,:);
    end
end

i_ids = 1;
data_label = {'experiment','fitting','sampling'}; %,'sample'};

% for i_data = 1:length(data.model_sim)
%     for i_data_type = 1:3
%         data_info.info_ligand{i_ids} = data.info_ligand{i_data};
%         data_info.info_dose_str{i_ids} = data.info_dose_str{i_data};
%         data_info.data_label{i_ids} = data_label{i_data_type};
%         i_ids = i_ids+1;
%     end
% end

data_sim = data;
metrics_sim = metrics_new;

load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

metrics_all = cell(1,length(metrics_sim)*3);

i_metric_index = 1;

metric_index_exp = [1:6,10:12,14:19]*2-1;
metric_index_fit = [1:6,10:12,14:19]*2;
data_index_exp = [1:6,10:12,14:19];
for i_metric_name = 1:length(metric_names)
    i_metric_index = 1;
    for i_metric = 1:length(metrics_sim)
        metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics{metric_index_exp(i_metric)}.(metric_names{i_metric_name})(:,:);
        data_info.info_ligand{i_metric_index} = data.info_ligand{data_index_exp(i_metric)};
        data_info.info_dose_str{i_metric_index} = data.info_dose_str{data_index_exp(i_metric)};
        data_info.data_label{i_metric_index} = 'experiment';
        i_metric_index = i_metric_index +1;
        
    end
    
    metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics_unstim{1}.(metric_names{i_metric_name})(:,:);
    data_info.info_ligand{i_metric_index} = 'none';
    data_info.info_dose_str{i_metric_index} = '0';
    data_info.data_label{i_metric_index} = 'experiment';
    
    i_metric_index = i_metric_index +1;
    
    for i_metric = 1:length(metrics_sim)
        metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics{metric_index_fit(i_metric)}.(metric_names{i_metric_name})(:,:);
        data_info.info_ligand{i_metric_index} = data.info_ligand{data_index_exp(i_metric)};
        data_info.info_dose_str{i_metric_index} = data.info_dose_str{data_index_exp(i_metric)};
        data_info.data_label{i_metric_index} = 'fitting';
        i_metric_index = i_metric_index +1;
        
    end
    
    metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metics_sim_unstim{1}.(metric_names{i_metric_name})(:,:);
    data_info.info_ligand{i_metric_index} = 'none';
    data_info.info_dose_str{i_metric_index} = '0';
    data_info.data_label{i_metric_index} = 'fitting';
    
    i_metric_index = i_metric_index +1;
    
    
    
    for i_metric = 1:length(metrics_sim)
        metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metrics_sim{i_metric}.(metric_names{i_metric_name})(:,:);
        data_info.info_ligand{i_metric_index} = data.info_ligand{data_index_exp(i_metric)};
        data_info.info_dose_str{i_metric_index} = data.info_dose_str{data_index_exp(i_metric)};
        data_info.data_label{i_metric_index} = 'sampling';
        i_metric_index = i_metric_index +1;
        
    end
    
    metrics_all{i_metric_index}.(metric_names{i_metric_name}) = metics_sim_unstim{1}.(metric_names{i_metric_name})(:,:);
    data_info.info_ligand{i_metric_index} = 'none';
    data_info.info_dose_str{i_metric_index} = '0';
    data_info.data_label{i_metric_index} = 'sampling';
    
    i_metric_index = i_metric_index +1;
    
end


vis_data_field = {'experiment','fitting','sampling'}; %,'sample'};
[collect_feature_vects,metrics] = calculate_codon_from_metric2023(data_info,metrics_all); %,  parameter
collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

for i_time_pts = 1:97
    time_list{i_time_pts} =strcat('time_points_',num2str(i_time_pts));
end


load('mutual_info_cal_data_example.mat')


index_ade = [3:3:15,16];
index_fitting = [19:3:31,32];
index_sampling = [35:3:47,48];


ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC','none'};

data_name = {'Ade','Fitting','Sampling'};
index_vec = {index_ade,index_fitting,index_sampling};

for i_data_set = 1:length(data_name)
    index_data = index_vec{i_data_set};
    
    for i_ligand = 1:length(ligand_all)
        nfkb(i_ligand).sc_metrics = struct();
        
        for i_time_pts =1:length(time_list)
            nfkb(i_ligand).id = ligand_all{i_ligand};
            nfkb(i_ligand).ids = ligand_all;
            nfkb(i_ligand).sc_metrics.(time_list{i_time_pts}) = metrics_all{index_data(i_ligand)}.time_series(:,i_time_pts);
        end
    end
    
    save(strcat('mutual_info_format_traj_single_ligand_',data_name{i_data_set},'.mat'),'nfkb')
    
end


codon_list = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};

for i_data_set = 1:length(data_name)
    index_data = index_vec{i_data_set};
    
    for i_ligand = 1:length(ligand_all)
        nfkb(i_ligand).sc_metrics = struct();
        
        for i_codon =1:length(codon_list)
            nfkb(i_ligand).id = ligand_all{i_ligand};
            nfkb(i_ligand).ids = ligand_all;
            % nfkb(i_ligand).sc_metrics.(time_list{i_time_pts}) = metrics_all{index_data(i_ligand)}.time_series(:,i_time_pts);
            nfkb(i_ligand).sc_metrics.(codon_list{i_codon}) = collect_feature_vects.(codon_list{i_codon}){index_data(i_ligand)};

        end
    end
    
    save(strcat('mutual_info_format_codon_single_ligand_unstim_',data_name{i_data_set},'.mat'),'nfkb')
    
end