% goal: plot heatmap of supriya's data, ade's data, fitting data, sampling data
% in the order of total activitiy, freq, osc power, etc.
% need to set to 1 to draw save figs 
%% notes
% Feature_Name	specifier	Codeword
% duration	2	Duration
% time2HalfMaxPosIntegral	-1	EarlyVsLate
% oscpower	0	OscVsNonOsc
% max_value	0	PeakAmplitude
% pos_pk1_amp	0	PeakAmplitude
% max_pos_pk1_speed	0	Speed
% pos_pk1_time	-1	Speed
% derivatives	2	Speed
% max_pos_integral	0	TotalActivity


% TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
% Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
% CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
% LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
% PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15

% goal: supriya's data, ade's data, fitting data, and sampled data for
% benchamarking
% compare the prediction of smapled data with supriya's data

% TNF 10ng/mL Pam3CSK 100ng/mL supriya: 03 index: 3
% TNF 10ng/mL CpG 100nM supriya: 03 index: 5
% TNF 10ng/mL LPS 10ng/mL supriya: 04 index: 8
% TNF 10ng/mL PolyIC 100ug/mL supriya: 04 index: 9
% Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
% Pam3CSK 100ng/mL LPS 10ng/mL supriya: 23 index: 26
% Pam3CSK 100ng/mL PolyIC 100ug/mL supriya: 13 index: 14
% CpG 100nM LPS 10ng/mL supriya: 23 index: 25
% CpG 100nM PolyIC 100ug/mL supriya: 21 index: 18
% LPS 10ng/mL PolyIC 100ug/mL supriya: 21 index: 19


%% initalization 
% clear all
addpath('./lib/')
addpath('./src/')
addpath('./bin/')

fig_save_path = '../SubFigures2023/';

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';


%% plot heatmap of ade's data

if 0
    clear data_to_draw data metrics
    
    % Ade's data:
    
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    data_index = [3;18;11;6;16];
    
    vis_data_field = {'exp_ade'};
    
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}
    data_to_draw.(vis_data_field{1}) = data.exp_mode_filter_nan(data_index);
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);
    
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)*2-1}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end

%% plot heatmap of fitting
if 0
    clear data_to_draw data metrics
    
    % ade's data
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    data_index = [3;18;11;6;16];
    vis_data_field = {'fitting'};
    
    clear data_to_draw
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}
    data_to_draw.(vis_data_field{1}) = data.pred_mode_filter_nan(data_index);
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);
    
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)*2}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end

%% plot heatmap of supriya's data
if 0
    clear data_to_draw data metrics
    
    % Supriya's data
    load('Supriya_data_metrics.mat')
    data_index = [4;27;22;6;15];
    vis_data_field = {'exp_supriya'};
    
    
    clear data_to_draw
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}
    data_to_draw.(vis_data_field{1}) = data.exp(data_index);
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);
    
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end

%% plot heatmap of sampling
if 0
    clear data_to_draw data metrics
    
    % sampling  data
    data_filename = 'Sim5_codon_all5dose_metric.mat';
    load(strcat(data_save_file_path_1,data_filename))
    data_index = [13;12;11;14;15];
    vis_data_field = {'sampling'};
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}

    
    metric_fields_name = fieldnames(metrics{1});
    for i_metric = 1:length(metrics)
        
        for i_metric_fields = 1:length(metric_fields_name)
            metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
        end
    end
    metrics = metrics_sample;
    
    clear data_to_draw
    data_to_draw.(vis_data_field{1}) = data.model_sim(data_index);
    for i_data = 1:length(data_to_draw.(vis_data_field{1}))
        data_to_draw.(vis_data_field{1}){i_data} = data_to_draw.(vis_data_field{1}){i_data}(1:9:end,:);
    end
    
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);

    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end

%% plot heatmap of dual ligand sampling
if 0
    clear data_to_draw data metrics
    
    % sampling  data
    data_filename = 'Sim5_codon_all5dose_metric.mat';
    load(strcat(data_save_file_path_1,data_filename))
    data_index = [1;2;3;4;10;9;6;8;5;7];

    vis_data_field = {'sampling_dual_ligand'};
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}

    
    metric_fields_name = fieldnames(metrics{1});
    for i_metric = 1:length(metrics)
        
        for i_metric_fields = 1:length(metric_fields_name)
            metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
        end
    end
    metrics = metrics_sample;
    
    clear data_to_draw
    data_to_draw.(vis_data_field{1}) = data.model_sim(data_index);
    for i_data = 1:length(data_to_draw.(vis_data_field{1}))
        data_to_draw.(vis_data_field{1}){i_data} = data_to_draw.(vis_data_field{1}){i_data}(1:9:end,:);
    end
    
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);

    for i_info_ligand = 1:length(data_to_draw.info_ligand)
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'polyIC','PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');  
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');      
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');      
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_PolyIC','LPS_PolyIC');
     end
    
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
% sampled data

end

%% plot heatmap of dual ligand supriya
if 0
    clear data_to_draw data metrics
    
    % Supriya's data
    load('Supriya_data_metrics.mat')
    data_index = [3;5;8;9;26;14;25;18;19];

    vis_data_field = {'exp_supriya'};
    
    
    clear data_to_draw
    order_name_vec = {'max_pos_integral','oscpower','oscfreq','peakfreq'};%,'pos_pk1_amp'}
    data_to_draw.(vis_data_field{1}) = data.exp(data_index);
    data_to_draw.info_ligand = data.info_ligand(data_index);
    data_to_draw.info_dose_str = data.info_dose_str(data_index);
    
        for i_info_ligand = 1:length(data_to_draw.info_ligand)
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'polyIC','PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');  
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');      
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');      
        data_to_draw.info_ligand{i_info_ligand} = replace(data_to_draw.info_ligand{i_info_ligand},'CpG_PolyIC','LPS_PolyIC');
        end
     
    for i_order_name = 1:length(order_name_vec)
        for i_data = 1:length(data_to_draw.(vis_data_field{1}))
            
            [~,data_to_draw.order{i_data}] = sort(metrics{data_index(i_data)}.(order_name_vec{i_order_name}),'descend');
        end
        filter_TNF = 0;
        plot_traj_heatmap_2023_05(data_to_draw,vis_data_field,fig_save_path,filter_TNF,order_name_vec{i_order_name})
    end
    
end



%% might be useful
%
%
% metrics = [metrics_supriya(supriya_index),... Supriya's single liagand data
%     metrics_ade(ade_index*2-1),... Ade's single liagand data
%     metrics_ade(ade_index*2),... Fitting Ade single liagand data
%     metrics_sample(sample_index),... Sampling single liagand data
%     metrics_supriya(supriya_dual_ligand_index),... Supriya's dual ligand data
%     metrics_sample(sample_dual_ligand_index)]; % sampling dual ligand data
%
% data_info.info_ligand = [data_supriya.info_ligand(supriya_index),... Supriya's single liagand data
%     reshape(data_ade.info_ligand(ade_index),1,[]),... Ade's single liagand data
%     reshape(data_ade.info_ligand(ade_index),1,[]),... Fitting Ade single liagand data
%     data_sample.info_ligand(sample_index),... Sampling single liagand data
%     data_supriya.info_ligand(supriya_dual_ligand_index),... Supriya's dual ligand data
%     data_sample.info_ligand(sample_dual_ligand_index)]; % sampling dual ligand data
%
% for i_info_ligand = 1:length(data_info.info_ligand)
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'polyIC','PolyIC');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');
%     data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_PolyIC','LPS_PolyIC');
% end
%
% data_info.info_dose_str = [data_supriya.info_dose_str(supriya_index),... Supriya's single liagand data
%     reshape(data_ade.info_dose_str(ade_index),1,[]),... Ade's single liagand data
%     reshape(data_ade.info_dose_str(ade_index),1,[]),... Fitting Ade single liagand data
%     data_sample.info_dose_str(sample_index),... Sampling single liagand data
%     data_supriya.info_dose_str(supriya_dual_ligand_index),... Supriya's dual ligand data
%     data_sample.info_dose_str(sample_dual_ligand_index)]; % sampling dual ligand data
%
% data_info.data_label(1:length(supriya_index)) = {'Exp_Supriya'};
% data_info.data_label(end+1:end+length(ade_index)) = {'Exp_Ade'};
% data_info.data_label(end+1:end+length(ade_index)) = {'Fitting'};
% data_info.data_label(end+1:end+length(sample_index)) = {'Sampling'};
% data_info.data_label(end+1:end+length(supriya_dual_ligand_index)) = {'Exp_Supriya_dual_ligand'};
% data_info.data_label(end+1:end+length(sample_dual_ligand_index)) = {'Sampling_dual_ligand'};
%
% [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023(data_info,metrics); %,  parameter
% metrics = metrics_new;
% %
% collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
% collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
% collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
% %
% save(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'),'collect_feature_vects','metrics','data_info')
% % goal: visualize supriya, ade, fitting, sampling (2 seprate plots) as
% % benchmarking; then visualize supriya dual ligand vs sampling.
