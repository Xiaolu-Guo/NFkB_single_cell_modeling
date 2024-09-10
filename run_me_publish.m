        % run_me_2023.m for NFkB parameter estimation
% addpath with subfolders: MACKtrack

%All_ligand_codon_2023_t33_cv_filtered_TNF.mat
% data.parameters_mode_nan
% data.exp_mode_filter_nan
% data.pred_mode_filter_nan
% are the TNF cv filtered trajectories

%% to delete:

% For Guo et al. Figure 5A, extended doses of each ligand simulated heatmaps
%
% tested 05/12/2024, Matlab 2020a

%%
clear all
vers = '202301';

run_sim_example = 0;
rescale_exp_data = 0;
rescale_supriya_exp_data_benchmark = 0;
cal_codon_diff_supriya_exp_data_benchmark = 0;
run_save_all_data_codon = 0;

draw_fig_traj_heatmap = 0;
draw_fig_para_distr = 0;
run_save_codon = 0;
draw_codon = 0;
draw_spearman_corr = 0;
draw_heatmap_para_scatter = 0;
save_para_codon_corr_python = 0;
visualize_para_codon_importance_python = 0;

transfer_codon_data_mutual_info_cal_format = 0;
plot_codon_cc_LPS_TNF = 0;
transfer_codon_data_machine_learning_format = 0;
visualize_stim_classification = 0;

run_sampling = 0;
run_sampling_fitting = 0;
draw_sampling = 0;
draw_sampling_codon = 0;
run_unstim_metric = 0;

%_supriya_data_sampling
draw_co_sti_sampling = 0;
draw_co_sti_CpG_polyIC_competetation = 0;

run_multi_sti = 0;
run_5_sti = 0;
run_scramble_dual = 0;
run_5_single_ligand = 0;

% intrinsic
draw_intrinsic_noise = 0;

%% initializing
debug_or_not = 0;

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
monolix_data_save_file_path = '../raw_data2023/SAEM_proj_2023/';
fig_save_path = '../SubFigures2023/';
data_example_format_path = './example_data_format/';


addpath('./lib/')
addpath('./src/')
addpath('./bin/')
addpath(genpath('./refs_codes/'));


if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

%% [tested] run & draw parameter sens analysis
% Figure S1, 0507 version: tested 05/09/2024
if 0
    ParameterScan_Sens_traj_04172024()
end

%% [tested] draw Heatmaps of traj, W-dist of signaling codon distri.
% Figure 2, 0507 version: tested 05/10/2024
if 0
    
    % Figure 2A
    draw_traj_heatmap_2023_03(monolix_data_save_file_path,fig_save_path)
    
    % Figure 2B-C
    draw_all_ligand_codon_distrib_202404(data_save_file_path,fig_save_path)
    
end

%% [tested] draw stimulation classification
% Figure 2D & Figure S4D : 0507 version: : tested 05/30/2024
if 0
    visualize_stim_classification_machine_learning_results
    ML_visual_sens_cross_data
end

%% [tested] draw metrics of good fitting, RMSD, signaling codon distri, and W-dist
% Figure S2: 0507 version: : tested 05/09/2024
if 0
    
    % Figure S2A RMSD distribution tested 06/18/2024
    draw_traj_RMSD_distrib_2024_06(fig_save_path)
    
    % Figure S2A
    draw_traj_RMSD_2023_05(fig_save_path)
    
    % Figure S2B
    draw_traj_RMSD_2024_05(fig_save_path)
    
    % Figure S2C-E: tested 05/10/2024
    draw_all_ligand_codon_distrib_202309(data_save_file_path,fig_save_path)
    
end

%% [to delete] [tested] save the parameter signaling codon for python codes: linear regression, RandomForest, XGBoost
% For Figure 3:
% save the data of parameters (X), and singnaling codon (y), and random
% variables as negative control. these dataset will be used to train the
% regression models in python. the outputs of the python codes have been
% saved in raw_data2023/singaling_codon_para/, and are visualized in the
% next two sections
if 0
    py_data_save_path = strcat(data_save_file_path,'singaling_codon_para_for_python/');
    
    save_csv_para_codon_2023_03(data_save_file_path,py_data_save_path,monolix_data_save_file_path)
end

%% [to delete] [tested] draw R-squared for regression models  and spearman corr and feature importance between para and codon
if 0 % Figure 3
    % Figure 3A 3C: 0507 version: : tested 05/12/2024
    visualize_para_codon_feature_importance_2024_04(data_save_file_path,py_data_save_path,fig_save_path)
    
    % Fgiure 3B: 0507 version: : tested 05/12/2024
    draw_spearman_corr_each_codon_2024_04(data_save_file_path,fig_save_path)
    
end

%% [tested] draw codon vs parameter: heatmap, spearman corrlation, feature importance
if 0 % Figure S3
    py_data_save_path = strcat(data_save_file_path,'singaling_codon_para/');
    
    % Figure S3A: 0507 version: : tested 05/12/2024
    para_cells_scatter_2023_06();
    
    % Figure S3B spearman correlation plot: 0507 version: : tested 05/12/2024
    draw_spearman_corr_each_codon_2023_06(data_save_file_path,fig_save_path)
    
    % Figure S3C: 0507 version: : tested 05/12/2024
    visualize_para_codon_feature_importance_2023_03(data_save_file_path,py_data_save_path,fig_save_path)
    
end

%% codon prep for MI calculation & machine learning format
% run till here
if 1% transfer_codon_data_mutual_info_cal_format % transfer the codon data to the mutula information calculation format
    
    if 0 % tested: 05/15/2024. tansfer sampling data
        data_save_file_path_1 = '../raw_data2023/';
        MI_file_save_path = '../raw_data2023/MI_single_ligand/';
        transfer_data_mutual_info_cal_format(data_save_file_path_1,MI_file_save_path)
    end
    

    if 0 % transfer 20doses
        % 'Sim3_codon_r5_metric.mat'; Pam3CSK
        % 'Sim3_codon_r4_metric.mat'; PolyIC
        % 'Sim3_codon_r3_metric.mat'; CpG
        % 'Sim3_codon_r2_metric.mat'; LPS
        % 'Sim3_codon_r1_metric.mat'; TNF
        
        data_filename = 'Sim3_codon_r5_metric.mat';
        save_sampling_20doses_for_channel_capacity(data_filename)
        scatter_codon_sampling_sim_cc
    end
    
    if 0 % transfer experiment data
        transfer_exp_data_mutual_info_cal_format
    end
end


if plot_codon_cc_LPS_TNF %debug for figure 5
    scatter_codon_exp_sampling_sim_cc
end

if transfer_codon_data_machine_learning_format % transfer the codon data to the mutula information calculation format
    transfer_data_machine_learning_cal_format
end

if 0 % transfer the doses codon data to the mutula information calculation format
    transfer_doses_codon_ML_cal_format
end

%% [to delete] [tested] run sampling dual ligand stim
% For Figure 4:
% run single-ligand and dual-ligand stimulation sampling
% calculate the corresponding signaling codons
% save as 'Sim5_codon_all5dose_metric.mat'
if 0
    Sim_dual_ligand_stim
end

%% [to delete] [tested] draw dual ligand sampling & exp traj heatmap, signaling codon distribution
% check  codon
if 0
    
    % Figure 4B: 0507 version: : tested 05/12/2024
    draw_dual_ligand_traj_heatmap_2024_05
    
    % Figure 4C: 0507 version: : tested 05/12/2024
    compare_sample_dual_ligand_codom_202405()
    
end

%% [to delete] [tested] run co-stim CpG-polyIC competetation
% For Figure 4F-G:
% run CpG-polyIC dual-ligand stimulation sampling with competition
% calculate the corresponding signaling codons
% save as 'Sim7_CpG_polyIC_competation_metric.mat'
if 0
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    % run all conditions
    % Initializing_sampling_supriya_data_2023
    
    Initializing_sampling_supriya_data_2023_05
    index_CpG_polyIC = [14];%1;2;6;7;11;
    input_paras.proj_num_vec = input_paras.proj_num_vec(index_CpG_polyIC);
    input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_CpG_polyIC);
    input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_CpG_polyIC);
    input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_CpG_polyIC);
    input_paras.Num_sample = 300;
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    save_metric_name = 'Sim7_CpG_polyIC_competation_metric.mat';
    cal_codon =1;
    
    para_sample_CpG_polyIC_Compete_202405(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    
end

%% [to delete] [tested] draw co-sti CpG-polyIC competetation
if  0
    data_file = 'Sim7_CpG_polyIC_competation_metric.mat';
    
    % Figure 4F: 0507 version: : tested 05/12/2024
    draw_traj_heatmap_co_sti_CpG_polyIC_competetation(data_file,fig_save_path,data_save_file_path) ;
    
    % Figure 4G: 0507 version: : tested 05/12/2024
    CpG_polyIC_compete_compare_supriya_202309(data_file)
end

%% [tested] draw single ligand sampling
% Figure S4: 0507 version: : tested 05/10/2024
if 0 % heatmap
    
    % Figure S4A: 0507 version: : tested 05/10/2024
    draw_traj_heatmap_diff_order_2023_08
    
    % Figure S4B-C: 0507 version: : tested 05/12/2024
    draw_all_ligand_sampling_codon_distrib_202306(data_save_file_path,fig_save_path)
    
end

%% [tested] [to check results] run extended doses for denoise
if 1
    if 0% run data
    cal_ext_dose_denoise
    end
    
    if 0
        MI_file_save_path = '../raw_data2023/';
        % transfer_ED_MI_format(data_save_file_path,MI_file_save_path)
        transfer_ED_MI_format_LPS(data_save_file_path,MI_file_save_path)
        
    end
    
    if 0 % LPS 5 dose cal
                MI_file_save_path = '../raw_data2023/';

        transfer_ED_MI_format_LPS_5doses(data_save_file_path,MI_file_save_path)
    end
    
end

%% [tested] draw for extended doses study of each ligand
% Figure 5A: 0507 version: : tested 05/10/2024
if 0
    ligand_fig_save_path = strcat(fig_save_path,'doses_20/');
    for i_r = 1:5
        
        data_filename = strcat('Sim3_codon_r',num2str(i_r),'_metric.mat');
        draw_sampling_traj_heatmap_2023_05(ligand_fig_save_path,data_filename)
    end
end

if 0 % [tested] 06/12/2024 draw 20 doses signaling codon distributions
    % 'Sim3_codon_r5_metric.mat'; Pam3CSK
    % 'Sim3_codon_r4_metric.mat'; PolyIC
    % 'Sim3_codon_r3_metric.mat'; CpG
    % 'Sim3_codon_r2_metric.mat'; LPS
    % 'Sim3_codon_r1_metric.mat'; TNF
    
    draw_sampling_20doses_codons_distrib()
       
end
    
%% run sampling, not sure what this is ??
if run_sampling
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    %    run_sampling_simulation(monolix_data_save_file_path,fig_save_path)
    %    TNFo_dual_para_sample(vers,data_save_file_path,Num_sample)
    %     Num_sample = 500;
    
    %   para_sample_sim(data_save_file_path,Num_sample)
    %     var_fold = [2,4,8];
    %     ligand_str_ = 'TNF';
    % run all conditions
    %     Initializing_sampling_supriya_data_2023
    input_paras.Num_sample = 10;
    % input_paras.proj_num_vec ={[6,3,5]
    %     [6,3,7]};
    % input_paras.proj_ligand_vec = {{'polyIC','TNF','CpG'}
    %     {'polyIC','TNF','Pam3CSK'}};
    % input_paras.proj_dose_str_vec = { {'100ug','10ng','100nM'}
    %     {'100ug','10ng','100ng'}};
    % input_paras.proj_dose_val_vec = {{100000,10,100}
    %     {100000,10,100}};
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    input_paras.proj_num_vec = {[2];[2];[2];
        [3];[3];[3];
        [4];[4];[4];
        [5];[5];[5];
        [6];[6];[6]};
    input_paras.Num_sample = [30;30;30;
        10;10;10;
        10;10;10;
        10;10;10;
        10;10;10];
    
    input_paras.proj_ligand_vec = {{'TNF'};{'TNF'};{'TNF'};
        {'LPS'};{'LPS'};{'LPS'};
        {'CpG'};{'CpG'};{'CpG'};
        {'PolyIC'}; {'PolyIC'}; {'PolyIC'};
        {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'}};
    input_paras.proj_dose_str_vec = {{'100pg/mL'};{'1ng/mL'};{'10ng/mL'};
        {'1ng/mL'};{'3ng/mL'};{'10ng/mL'};
        {'33nM'};{'100nM'};{'333nM'};
        {'10ug/mL'};{'33ug/mL'};{'100ug/mL'};
        {'10ng/mL'};{'100ng/mL'};{'1ug/mL'}};
    input_paras.proj_dose_val_vec = {{0.1};{1};{10};
        {1};{3};{10};
        {33};{100};{333};
        {10000};{33000};{100000};
        {10};{100};{1000}};
    
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    cal_codon =1;
    save_metric_name = strcat('ligand5_sim_alldose_cellnumber10_codon_metric.mat');
    
    para_sample_sim_sti_doses_var_2023_03(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)
    
    %     for i_var_fold = 1:length(var_fold)
    %         para_sample_sim_sti_doses_small_var(data_save_file_path,Num_sample,ligand_str_,var_fold(i_var_fold))
    %     end
    %
    %     vers_fig = 'v1';
    %
    %     single_or_dual = 'dual';
    %     draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path)
    
end

%% sampling fitting run, not sure what this is for ??
if run_sampling_fitting
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    % draw_traj_heatmap_2023_03(monolix_data_save_file_path,fig_save_path)
    
    
    %     input_paras.proj_num_vec = {[2];[2];[2];
    %         [3];[3];[3];
    %         [4];[4];[4];
    %         [5];[5];[5];
    %         [6];[6];[6]};
    %     input_paras.Num_sample = [30;30;30;
    %         10;10;10;
    %         10;10;10;
    %         10;10;10;
    %         10;10;10];
    %     input_paras.proj_ligand_vec = {{'TNF'};{'TNF'};{'TNF'};
    %         {'LPS'};{'LPS'};{'LPS'};
    %         {'CpG'};{'CpG'};{'CpG'};
    %         {'PolyIC'}; {'PolyIC'}; {'PolyIC'};
    %         {'Pam3CSK'};{'Pam3CSK'};{'Pam3CSK'}};
    %     input_paras.proj_dose_str_vec = {{'100pg/mL'};{'1ng/mL'};{'10ng/mL'};
    %         {'1ng/mL'};{'3ng/mL'};{'10ng/mL'};
    %         {'33nM'};{'100nM'};{'333nM'};
    %         {'10ug/mL'};{'33ug/mL'};{'100ug/mL'};
    %         {'10ng/mL'};{'100ng/mL'};{'1ug/mL'}};
    %     input_paras.proj_dose_val_vec = {{0.1};{1};{10};
    %         {1};{3};{10};
    %         {33};{100};{333};
    %         {10000};{33000};{100000};
    %         {10};{100};{1000}};
    
    input_paras.proj_num_vec = {[2];[2];[2]};
    input_paras.Num_sample = [30;30;30];
    input_paras.proj_ligand_vec = {{'TNF'};{'TNF'};{'TNF'}};
    input_paras.proj_dose_str_vec = {{'100pg/mL'};{'1ng/mL'};{'10ng/mL'}};
    input_paras.proj_dose_val_vec = {{0.1};{1};{10}};
    
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    cal_codon =1;
    save_metric_name = strcat('ligand5_sim_alldose_cellnumber10_codon_metric.mat');
    
    % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)
    
    para_sample_fitting_alldose_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras,cal_codon,save_metric_name)
end

%
% TNF_color = [119 180 202]/255;
% CpG_color = [137 180 66]/255;
% P3CSK_color = [229 129 56]/255;
% LPS_color = [222 78 66]/255;
% PolyIC_cclor = [101 77 123]/255;

%% [tested] rescaling the data to SI
if 0 % rescale_exp_data and save SAEM format
    % tested 05/15/2024
    data_main
    % cd(pwd)
end

%% [tested] unstim codon cal
if 0
    % run simulation for unstimulated case: : tested 05/15/2024
    % to get Sim_unstim_fitting_alldose_r2_codon_metric.mat
    Sim_unstim_202405
    % save the data format into MI calculation format: tested 05/15/2024
    % to get
    % mutual_info_format_codon_single_ligand_unstim_Sampling_20230614.mat
    save_unstim_metric_for_channel_capacity_2024
end

%% [check later] rescaling the data for calculating benchmark for fitting
if 0 % rescale_supriya_exp_data_benchmark
    currentFolder = pwd;
    addpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_data/')
    data_main_benchmark_supriya_data
    % cd(pwd)
end

if 0 % cal_codon_diff_supriya_exp_data_benchmark
    cal_codon_benchmark_supriya_data
end

%% [tested] single cell specificity
if 1
    
    %% [tested] [run data] 5 single ligand stim: matching
    if 0  % tested 05/15/2024, Matlab 2020a
        %run_5_single_ligand: to get 'Sim8_5_signle_ligand_codon_metric_r3.mat'
        % doses info:
        % ligands: {{'TNF'};
        %     {'LPS'};
        %     {'CpG'};
        %     {'PolyIC'};
        %     {'Pam3CSK'}};
        
        % dose_str: {{'10ng/mL'};
        %     {'10ng/mL'};
        %     {'100nM'};
        %     { '100ug/mL'};
        %     {'100ng/mL'}};
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample = 200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        %save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r2.mat';
        for i_r = 3%:7
            save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_sim_single_ligand_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
    end
    
    %% [tested] [run data] single cell specificity for IkBa-/- match
    if 0 % tested 05/15/2024, Matlab 2020a
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample = 200;%5 ligand 200
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        
        cal_codon =1;
        
        if 1 % 03/22/2024 params6 = 0.25 * params6_wt;
            for  i_r =1
                save_metric_name = strcat('Sim16_IkBao_matching_5_signle_ligand_codon_metric_p25x_r',num2str(i_r),'.mat');
                % params6 = 0.1 * params6_wt;
                para_sample_fitting_match_sim_SRS_IkBao_202401(0.25, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
            end
        end
        
        if 0% params6 = 0.1 * params6_wt;
            for  i_r =1
                save_metric_name = strcat('Sim16_IkBao_matching_5_signle_ligand_codon_metric_p1x_r',num2str(i_r),'.mat');
                % params6 = 0.1 * params6_wt;
                para_sample_fitting_match_sim_SRS_IkBao_202401(0.1, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
            end
        end
        
        i_r =1;
        
        % params6 = 0.01 * params6_wt;
        if 0
            save_metric_name = strcat('Sim16_IkBao_matching_5_signle_ligand_codon_metric_p01x_r',num2str(i_r),'.mat');
            para_sample_fitting_match_sim_SRS_IkBao_202401(0.01, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        
        % params6 = 0 * params6_wt;
        if 0
            save_metric_name = strcat('Sim16_IkBao_matching_5_signle_ligand_codon_metric_0x_r',num2str(i_r),'.mat');
            para_sample_fitting_match_sim_SRS_IkBao_202401(0, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
    end
    
    %% [tested] [draw figure] draw signle cell confusion
    if 1 % Fgiure 7
        
        % figure 5
        scSRS_corr_heterogeneity
        
        % figure 5
        if 0
        draw_single_cell_confusion_20240404
        end
        % figure S5
        if 0
        Figure_S7A_single_cell_specifcity
        end
    end
    
end

%% [tested] Application: MI within the network
% run simulation, transfer into MI calculation format
% tested 05/15/2024
if 0
    
    data_save_file_path = '../raw_data2023/simulation_denoise/';
    MI_file_save_path = '../raw_data2023/MI_denoise/';
    
    if 0 % [tested] single ligand sampling for Sjroen syndrome
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        %Supriya's data single ligand dose
        %     input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        %     input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample =1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        cal_codon =1;
        
        if 1
            
            for i_r =1:2
                
                save_metric_name = strcat('Sim5_SS_codon_metric_p25x_r',num2str(i_r),'.mat');
                para_fitting_sim_IkBao_202401(0.25, data_save_file_path,input_paras,cal_codon,save_metric_name);
                
            end
        end
        
        if 0 % not run yet, 04/01/2024, there is some results in the denoise analysis, that might be checked
            for i_r =1:2
                
                save_metric_name = strcat('Sim5_WT_codon_metric_1x_r',num2str(i_r),'.mat');
                para_fitting_sim_IkBao_202401(1, data_save_file_path,input_paras,cal_codon,save_metric_name);
                
            end
        end
        
        if 0
            for i_r =1:2
                
                save_metric_name = strcat('Sim5_SS_codon_metric_p5x_r',num2str(i_r),'.mat');
                para_fitting_sim_IkBao_202401(0.5, data_save_file_path,input_paras,cal_codon,save_metric_name);
                
            end
        end
        
    end
    
    if 1% [tested] 06/18 cal para CV
        
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        data_label = {'TNF-syn','TNF-deg',...
            'LPS-syn','LPS-deg','LPS-endo',...
            'CpG-syn','CpG-deg','CpG-endo',...
            'PolyIC-syn','PolyIC-deg','PolyIC-endo',...
            'Pam3CSK-syn','Pam3CSK-deg','adp','core-timed1','core-timed2','core-NFkBtot'};
       
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {};%core IKK-NFkB-IKBa
        var_input.paravals = []';
        
        var_input.specienames = {};
        var_input.specievals = [];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Simxx_sampling_cal_CV','.mat');
        rcp_cvs = [];
        adp_cvs = [];
        core_cvs = [];
        for i=1:10
            % [rcp_cvs_rpc,adp_cvs_rpc,core_cvs_rpc] = para_sampling_cal_CV(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
            [rcp_cvs_rpc,adp_cvs_rpc,core_cvs_rpc] = para_sampling_cal_CV_barplot(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);

            rcp_cvs = [rcp_cvs,rcp_cvs_rpc];
            adp_cvs =[adp_cvs,adp_cvs_rpc];
            core_cvs = [core_cvs,core_cvs_rpc];
        end
        
        if 1 % bar plot of paramter cv
            figure(1)
            paperpos=[0,0,300,100]*1.5;
            papersize=[300 100]*1.5;
            draw_pos=[10,10,290,90]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            core_cvs_new = [core_cvs(1:3,:),core_cvs(4:6,:),core_cvs(7:9,:),core_cvs(10:12,:),core_cvs(13:15,:)];
            % cvs_all = [rcp_cvs;adp_cvs;core_cvs];
            % rcp_cvs = rcp_cvs(:,[
            bar_means_all = mean(rcp_cvs, 2);
            bar_stddevs_all = std(rcp_cvs, 0, 2);
            bar_means_all(end+1) = mean(adp_cvs(:));
            bar_stddevs_all(end+1) = std(adp_cvs(:));
            bar_means_all(end+1:end+3) = mean(core_cvs_new,2);
            bar_stddevs_all(end+1:end+3) = std(core_cvs_new, 0, 2);
            
            mean(mean(rcp_cvs, 2))
            mean(mean(adp_cvs(:)))
            mean(mean(core_cvs_new,2))
%             data_label = {'TNF-syn','TNF-deg',...
%             'LPS-syn','LPS-deg','LPS-endo',... 3-5
%             'CpG-syn','CpG-deg','CpG-endo',... 6-8
%             'PolyIC-syn','PolyIC-deg','PolyIC-endo',... 9-11
%             'Pam3CSK-syn','Pam3CSK-deg','adp','core-timed1','core-timed2','core-NFkBtot'};% 12, 13, 14 adp, 15-17
            data_label_index = [1,2,12,13,6,7,8,3,4,5,9,10,11,15,16,17];
            
            bar_means = bar_means_all(data_label_index);
            bar_stddevs = bar_stddevs_all(data_label_index);
            data_label = data_label(data_label_index);
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize)
            c = categorical(data_label);
            c = reordercats(c,data_label);
            bar(c,bar_means,0.4,'EdgeColor',[0 0 0],'LineWidth',0.5);hold on; % Adjust bar width as needed
            
            numGroups = length(bar_means);
            numBars = length(bar_means);
            groupWidth = min(0.8, numBars/(numBars + 1.5));
            x = (1:numGroups) - groupWidth/2 + (2*numBars-1) * groupWidth / (2*numBars); % Adjust the position
            
            % Add error bars
            errorbar(x, bar_means, bar_stddevs, 'k', 'linestyle', 'none');
            
            ax2 = gca;
            % ytickformat(ax2, '%g%%');
            % ylim([75,100])
%             xticklabels({});
%             yticklabels({});
            % title(rmsd_cas)
            set(gca,'fontsize',7,'XColor','k','YColor','k','FontName','Arial')
            
            
            saveas(gcf,strcat(fig_save_path,'cvs_bar_rcp_adp_core'),'epsc');
            
            close
            
        end
        
        
        if 0 % violin plot of distribution
            figure(1)
            paperpos=[0,0,100,100]*1.5;
            papersize=[100 100]*1.5;
            draw_pos=[10,10,90,90]*1.5;
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            
            y =  {rcp_cvs;adp_cvs;core_cvs};
            z =  {rcp_cvs;adp_cvs;core_cvs};
            % subplot(1,length(vis_data_field),i_data_field)
            
            al_goodplot_pair_RMSD_diff_size(y,[],0.5,ones(length(y),1)*[0 0 255]/255 ,'left',[],std(cell2mat(y))/2500);
            al_goodplot_pair_RMSD_diff_size(z,[],0.5,ones(length(z),1)*[0 0 255]/255,'right',[],std(cell2mat(z))/2500);
            
            xlim([0.4 3.6])
            
            xticks([1:3])
            xticklabels({})
            %title({strcat('K_{d,NFkB} =',num2str(params.Kd),', K_{d,p38} =',num2str(params.Kdp38))})
            
            ylim([0,2]);
            %     for i_x = 1:15
            %         plot([i_x,i_x],[0,5],'--','Color','k');hold on
            %     end
            set(gca,'fontsize',14,'fontname','Arial');
            %%%% saveas(gcf,strcat(fig_save_path,'PairRMSD_distrib_exp_',vers_savefig),'epsc');
            
            saveas(gcf,strcat(fig_save_path,'cvs_distrib_rcp_adp_core'),'epsc');
            
            close
            
        end
    


    end
    
    
    if 0 % [tested] no noise
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params52n2','params65n2',...
            'params99','params101'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            1,1,...
            0.4,0.4]';
        
        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim24_no_heterogeneity_codon_metric',num2str(i_r),'.mat');
        
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
    end
    
    if 0 %[tested] only TAK noise Sim22
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params99','params101'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            0.4,0.4]';
        
        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim22_TAK1_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
    end
    
    if 0 % [tested] only core module noise
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...LPS
            'params52n2','params65n2'};%core IKK-NFkB-IKBa
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            1,1]';
        
        var_input.specienames = {};
        var_input.specievals = [];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim23_NFkB_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
    end
    
    if 0 % [tested] only rcp noise
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params52n2','params65n2','params99','params101'};%LPS
        var_input.paravals = [1,1,0.4,0.4]';
        
        var_input.specienames = {'NFkB'};
        var_input.specievals = [0.08];
        
        cal_codon =1;
        
        for i_r =1:2
            save_metric_name = strcat('Sim21_reduce_TAK1ac_NFkB_heterogeneity_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end
        
        
    end
    
    if 0 % [tested] TAK reduce heterogeneity
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params52n2','params65n2'};%LPS
        var_input.paravals = [1,1]';
        
        var_input.specienames = {};
        var_input.specievals = [];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim20_reduce_TAK1ac_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
        
    end
    
    
    if 0 %% [tested] receptor reduce heterogeneity
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35'};%LPS
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05]';
        
        var_input.specienames = {};
        var_input.specievals = [];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim19_reduce_rcpt_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
        
    end
    
    if 0 %% [tested] receptor reduce heterogeneity for IkBas/s
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample =1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params93','params88','params85',...CpG
            'params83','params79','params77',...PolyIC
            'params75','params68',...Pam
            'params64','params61','params58','params54',...TNF
            'params44','params40','params36','params35',...%LPS
            'params6'};% NFkB regulated IkBa transcriptional rate
        var_input.paravals = [1.60E-03,0.015,2e-6,...
            0.0007,0.04,3e-6,...
            0.004,1e-6,...
            0.125,0.125,0.125,8.224e-06,...
            0.012,0.065681,0.065681,5.25E-05,...
            6e-05*0.25]';
        
        var_input.specienames = {};
        var_input.specievals = [];
        
        
        cal_codon =1;
        
        i_r =1;
        save_metric_name = strcat('Sim19_IkBas_reduce_rcpt_heterogeneity_codon_metric',num2str(i_r),'.mat');
        para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
        
    end
    
    if 0 % [tested] control: wt, IkB-/- % only use wt.
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 5;%1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        cal_codon =1;
        
        para_vec = [6e-05;
            6e-05*0.1;
            6e-05*0.01;
            6e-05*0];
        
        ncpu=2;
        pc=parcluster('local');
        pc.NumThreads=2;%
        parpool(pc,ncpu)
        var_input = cell(4);
        
        for i_r =1:size(para_vec,1)
            var_input{i_r}.paranames = {'params6'};
            var_input{i_r}.paravals = para_vec(i_r,:)';
            var_input{i_r}.specienames = {};
            var_input{i_r}.specievals = [];
            
            save_metric_name = strcat('Sim18_wt_IkBo_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input{i_r}, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end
        
        delete(gcp)
        
    end
    
    if 0 %[tested]  generate data, reducing NFkB var, using representaive, vs 9 random picked selected cells
        % purpose: assuming that all information are transimitted into NFkB
        
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF'},{ 'LPS'},{ 'CpG'},{ 'PolyIC'},{ 'Pam3CSK'}};
        % Supriya's data single ligand dose
        % input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'}};
        % input_paras.proj_dose_val_vec = {{10},{10},{100},{100000},{100}};
        
        %ade's data highest dose
        input_paras.proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'333nM'},{'100ug/mL'},{'1000ng/mL'}};
        input_paras.proj_dose_val_vec = {{10},{10},{333},{100000},{1000}};
        input_paras.Num_sample = 1000;%5 ligand
        input_paras.proj_num_vec ={[2],[3],[4],[5],[6]};
        
        var_input.paranames = {'params99','params101'};
        var_input.specienames = {'NFkB'};
        
        para_vec = [0.4,0.4;
            0.336004000000000,0.633038000000000;
            0.0673386000000000,0.0619458000000000;
            0.0157939000000000,0.0216838000000000;
            0.217826000000000,0.417812000000000;
            0.424314000000000,0.445427000000000;
            0.0115429000000000,0.0258170000000000;
            0.0802117000000000,0.0685842000000000;
            0.493319000000000,0.669808000000000;
            0.0241970000000000,0.0200540000000000];
        
        specie_vec = [0.08;
            0.0564181000000000;
            0.242900000000000;
            0.0701397000000000;
            0.0497623000000000;
            0.114479000000000;
            0.0974905000000000;
            0.0421724000000000;
            0.0721400000000000;
            0.0858496000000000];
        
        
        cal_codon =1;
        
        ncpu=2;
        pc=parcluster('local');
        pc.NumThreads=1;%
        parpool(pc,ncpu)
        var_input = cell(10);
        
        parfor i_r =1:size(para_vec,1)
            var_input{i_r}.paranames = {'params99','params101'};
            var_input{i_r}.specienames = {'NFkB'};
            var_input{i_r}.paravals = para_vec(i_r,:)';
            var_input{i_r}.specievals = specie_vec(i_r,:)';
            
            save_metric_name = strcat('Sim17_reduce_core_heterogeneity_codon_metric',num2str(i_r),'.mat');
            para_fitting_sim_alter_var_202402(var_input{i_r}, data_save_file_path,input_paras,cal_codon,save_metric_name);
        end
        
        delete(gcp)
        
    end
    
    if 0 % [tested]
        transfer_diff_var_codon_mutual_info_format_0331(data_save_file_path,MI_file_save_path)
    end
    
    
end

% before 01/31/2024
% transfer time trajectory into MI format
if 0
    transfer_sample_traj_into_mutual_info_format
end

if 0 % plot the trajs and codon distribution
    transfer_sample_traj_into_mutual_info_format
end

%% publish: not figures, numbers in the manuscript

%% [supplementary materials] benchmarking different methods combine receptor modules
if 0
    %% [tested] [supplementary materials] 5 single ligand stim: weighted matching
    if 0 % tested 05/15/2024, Matlab 2020a
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        save_metric_name = 'Sim15_5_signle_ligand_codon_metric.mat';
        
        % ncpu = 5;
        % pc = parcluster('local');
        % pc.NumThreads = 2;
        % parpool(pc,ncpu)
        % par
        for i_r = 2 :6
            
            save_metric_name = strcat('Sim15_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_weight_match_sim_single_ligand_var_2023_12(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        % delete(gcp)
    end
    
    %% [tested] [supplementary materials] 5 single ligand stim: scramble
    if 0 % tested 05/15/2024, Matlab 2020a
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        %    save_metric_name = 'Sim14_5_signle_ligand_codon_metric_r2.mat';
        
        % ncpu=5;
        % pc=parcluster('local');
        % pc.NumThreads=2;%
        % parpool(pc,ncpu)
        %
        % par
        for  i_r =3:7
            save_metric_name = strcat('Sim14_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_sim_single_ligand_scramble_2024(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        % delete(gcp)
        
    end
    
    %% [check later] [supplementary materials] single cell compare draw replicates
    if 0
        %draw_5_single_ligand_sampling_exp_codon_v202401
        draw_5_single_ligand_sampling_exp_codon_v20240120
    end
    
end

%% [to delete] calculate responder ratio: from supriya's data
if 0
    calculation_responder_5sigma_Supriyadata()
end

%% save all data codon: dealing with data
if run_save_all_data_codon
    save_all_data_codon
end

%% calculate signaling codon, and filter out all low variation of TNF oscillation
% data prepare for all codon visulization
% tested 05/10/2024
run_save_codon = 0;
if run_save_codon
    %run_all_ligand_codon_2023(monolix_data_save_file_path,data_save_file_path)
    run_all_ligand_codon_2023_03(monolix_data_save_file_path,data_save_file_path)
    % run_diff_ligand_codon(monolix_data_save_file_path,data_save_file_path)
end