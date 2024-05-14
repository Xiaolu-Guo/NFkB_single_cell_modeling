% run_me_2023.m for NFkB parameter estimation
% addpath with subfolders: MACKtrack

%All_ligand_codon_2023_t33_cv_filtered_TNF.mat
% data.parameters_mode_nan
% data.exp_mode_filter_nan
% data.pred_mode_filter_nan
% are the TNF cv filtered trajectories

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

%% [tested] run parameter sens analysis 
% Figure S1, 0507 version: tested 05/09/2024
if 0
    ParameterScan_Sens_traj_04172024()
end

%% [tested] Heatmaps of traj, W-dist of signaling codon distri.
% Figure 2, 0507 version: tested 05/10/2024
if 0
    
    % Figure 2A
    draw_traj_heatmap_2023_03(monolix_data_save_file_path,fig_save_path) 

    % Figure 2B-C
    draw_all_ligand_codon_distrib_202404(data_save_file_path,fig_save_path)

end

%% [tested] visualize stimulation classification
% Figure 2D & Figure S4D : 0507 version: : tested 05/11/2024
if 0
    visualize_stim_classification_machine_learning_results
end

%% [tested] metrics of good fitting, RMSD, signaling codon distri, and W-dist
% Figure S2: 0507 version: : tested 05/09/2024
if 0
     
    % Figure S2A
    draw_traj_RMSD_2023_05(fig_save_path) 
    
    % Figure S2B
    draw_traj_RMSD_2024_05(fig_save_path) 

    % Figure S2C-E: tested 05/10/2024
    draw_all_ligand_codon_distrib_202309(data_save_file_path,fig_save_path)
    
end

%% [tested] save the parameter signaling codon for python codes: linear regression, RandomForest, XGBoost
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

%% [tested] R-squared for regression models  and spearman corr and feature importance between para and codon
if 0 % Figure 3    
    % Figure 3A 3C: 0507 version: : tested 05/12/2024
    visualize_para_codon_feature_importance_2024_04(data_save_file_path,py_data_save_path,fig_save_path)

    % Fgiure 3B: 0507 version: : tested 05/12/2024
    draw_spearman_corr_each_codon_2024_04(data_save_file_path,fig_save_path)
     
end

%% [tested] codon vs parameter: heatmap, spearman corrlation, feature importance
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
if transfer_codon_data_mutual_info_cal_format % transfer the codon data to the mutula information calculation format
    if 0 % tansfer sampling data
        transfer_data_mutual_info_cal_format
    end
    if 0 % transfer 20doses
        % 'Sim3_codon_r5_metric.mat'; Pam3CSK
        % 'Sim3_codon_r4_metric.mat'; PolyIC
        % 'Sim3_codon_r3_metric.mat'; CpG
        % 'Sim3_codon_r2_metric.mat'; LPS
        % 'Sim3_codon_r1_metric.mat'; TNF
        
        data_filename = 'Sim3_codon_r5_metric.mat';
        % save_sampling_20doses_for_channel_capacity(data_filename)
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


%% [tested] run sampling dual ligand stim 
% For Figure 4: 
% run single-ligand and dual-ligand stimulation sampling
% calculate the corresponding signaling codons
% save as 'Sim5_codon_all5dose_metric.mat'
if 0
    Sim_dual_ligand_stim
end

%% [tested] draw dual ligand sampling & exp traj heatmap, signaling codon distribution
% check  codon
if 0
    
    % Figure 4B: 0507 version: : tested 05/12/2024
    draw_dual_ligand_traj_heatmap_2024_05
    
    % Figure 4C: 0507 version: : tested 05/12/2024
    compare_sample_dual_ligand_codom_202405()
    
end

%% [tested] run co-stim CpG-polyIC competetation 
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

%% [tested] draw co-sti CpG-polyIC competetation    
if  0
    data_file = 'Sim7_CpG_polyIC_competation_metric.mat';
    
    % Figure 4F: 0507 version: : tested 05/12/2024
    draw_traj_heatmap_co_sti_CpG_polyIC_competetation(data_file,fig_save_path,data_save_file_path) ;
    
    % Figure 4G: 0507 version: : tested 05/12/2024
    CpG_polyIC_compete_compare_supriya_202309(data_file)
end

%% [tested] single ligand sampling
% Figure S4: 0507 version: : tested 05/10/2024
if 0 % heatmap
    
    % Figure S4A: 0507 version: : tested 05/10/2024
    draw_traj_heatmap_diff_order_2023_08
    
    % Figure S4B-C: 0507 version: : tested 05/12/2024
    draw_all_ligand_sampling_codon_distrib_202306(data_save_file_path,fig_save_path)
    
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

%% rescaling the data to SI
if rescale_exp_data
    currentFolder = pwd;
    addpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_data/')
    data_main
    % cd(pwd)
end

%% unstim codon cal
if run_unstim_metric
    %     unstim_metric
    % updated version 06/14/2023
    save_unstim_metric_for_channel_capacity
end

%% rescaling the data for calculating benchmark for fitting
if rescale_supriya_exp_data_benchmark
    currentFolder = pwd;
    addpath('/Users/admin/Documents/my document/Postdoc projects/MatlabCodes/NFkB_data/')
    data_main_benchmark_supriya_data
    % cd(pwd)
end

if cal_codon_diff_supriya_exp_data_benchmark
    cal_codon_benchmark_supriya_data
end

%% multi stimulation
if run_multi_sti
    monolix_data_save_file_path = '../SAEM_proj_2023/';
    % run all conditions
    % Initializing_sampling_supriya_data_2023
    
    Initializing_sampling_supriya_data_2023_05
    input_paras.Num_sample = 2;%3 ligand
    % input_paras.proj_num_vec ={[6,3,5]
    %     [6,3,7]};
    % input_paras.proj_ligand_vec = {{'polyIC','TNF','CpG'}
    %     {'polyIC','TNF','Pam3CSK'}};
    % input_paras.proj_dose_str_vec = { {'100ug','10ng','100nM'}
    %     {'100ug','10ng','100ng'}};
    % input_paras.proj_dose_val_vec = {{100000,10,100}
    %     {100000,10,100}};
    % TNF 2, LPS 3, CpG 4, PolyIC 5, Pam 6
    index_multi_ligand = [21;22];%1;2;6;7;11;
    %    index_multi_ligand = [3];%1;2;6;7;11;
    
    input_paras.proj_num_vec = input_paras.proj_num_vec(index_multi_ligand);
    input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_multi_ligand);
    input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_multi_ligand);
    input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_multi_ligand);
    
    save_metric_name = 'Sim6_debug_3ligands_codon_metric.mat';
    cal_codon =1;
    para_sample_fitting_sim_sti_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    
    %     vers_fig = '2023';
    %     single_or_dual = 'dual';
    %     para_sample_sim(data_save_file_path,Num_sample)
    %     draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path)
end

%% all possible combinatorial ligands stim
run_all_comb_sti = 0;
if run_all_comb_sti
    
    if 0% 2-3-4-5 ligands
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        Initialize_combinatorial_ligands
        
        save_metric_name = 'Sim9_all_Comb_ligands_weighted_match_codon_metric.mat';
        cal_codon =1;
        % matched
        %para_sample_fitting_sim_sti_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        % weighted matche
        para_sample_fitting_weighted_match_sim_2024(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    end
    
    if 0 % single ligands
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        Initialize_combinatorial_single_ligands
        
        save_metric_name = 'Sim10_all_single_Comb_ligands_codon_metric.mat';
        cal_codon =1;
        para_sample_fitting_sim_sti_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    end
    
end

draw_all_comb_sti = 0;
if draw_all_comb_sti
    draw_combinatorial_ligand
end

%% single cell specificity
if 1
    % doses are based one supriya's dataset
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
    
    %% 5 single ligand stim: matching
    if 0 %run_5_single_ligand
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        %save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r2.mat';
        for i_r = 3:7
            save_metric_name = strcat('Sim8_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_sim_single_ligand_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
    end
    
    %% 5 single ligand stim: weighted matching
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        save_metric_name = 'Sim15_5_signle_ligand_codon_metric.mat';
        
        % ncpu = 5;
        % pc = parcluster('local');
        % pc.NumThreads = 2;
        % parpool(pc,ncpu)
        % par
        for i_r = 2:6
            
            save_metric_name = strcat('Sim15_5_signle_ligand_codon_metric_r',num2str(i_r),'.mat');
            cal_codon =1;
            para_sample_fitting_weight_match_sim_single_ligand_var_2023_12(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        end
        % delete(gcp)
    end
    
    %% 5 single ligand stim: scramble
    if 0
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
    
    %% 5 single ligand stim: sample from original distrib, how to match? seems the minimum dist as match? should randomly match?
    % to do: try random match
    if 0
        Sim2_sampling_2024
    end
    
    %% single cell compare draw replicates
    if 0
        %draw_5_single_ligand_sampling_exp_codon_v202401
        draw_5_single_ligand_sampling_exp_codon_v20240120
    end
    
    %% single cell specificity for IkBa-/- match
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand 200
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        
        cal_codon =1;
        
        if 0 % 03/22/2024 params6 = 0.25 * params6_wt;
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
    
    %% single cell specificity for IkBa-/- weighted match
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec ={{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample =200;%5 ligand
        input_paras.proj_num_vec ={[2,3,4,5,6]};
        
        cal_codon =1;
        
        if 0% params6 = 0.1 * params6_wt;
            for  i_r =1
                save_metric_name = strcat('Sim16_IkBao_5_signle_ligand_codon_metric_p1x_r',num2str(i_r),'.mat');
                % params6 = 0.1 * params6_wt;
                para_sample_fitting_weight_match_sim_SRS_IkBao_202401(0.1, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
            end
        end
        
        i_r =1;
        
        % params6 = 0.01 * params6_wt;
        save_metric_name = strcat('Sim16_IkBao_5_signle_ligand_codon_metric_p01x_r',num2str(i_r),'.mat');
        para_sample_fitting_weight_match_sim_SRS_IkBao_202401(0.01, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        
        % params6 = 0 * params6_wt;
        save_metric_name = strcat('Sim16_IkBao_5_signle_ligand_codon_metric_0x_r',num2str(i_r),'.mat');
        para_sample_fitting_weight_match_sim_SRS_IkBao_202401(0, data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        
    end
    
    %% draw signle cell confusion
    if 0
        draw_single_cell_confusion_20240404
        %draw_single_cell_confusion_20240319
    end
            
    %% draw single cell specificity heatmaps
    if 0 % step 2 Figure 7
        draw_single_cell_specifcity_matched_20240213
    end
    
    %% draw single cell specificity clusters for IkBo-/-
    if 0
        draw_single_cell_specifcity_matched
    end
    
    %% transfer all the codon into MI cal format
    if 0
        transfer_sample_SRS_data_MI_format
    end
    
    %% tranfer into para vs low-med-high sepcificity categories mat, for python machine learning
    if 0
        addpath('./python/')
        py_data_save_path = './python/';
        save_csv_para_srs_2024(data_save_file_path,py_data_save_path)
    end
    
    %% plot para distribution by cluster number
    if 0
        draw_paraDistrib_by_clusterNo_2024(data_save_file_path)
    end
    
    %% from here to be deleted
    % dual ligand match, single ligand reponse to be deleted??
    run_2_single_ligand = 0;
    
    if run_2_single_ligand
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        Initialize_combinatorial_ligands
        input_paras_all = input_paras;
        clear input_paras
        proj_num_vec_all = input_paras_all.proj_num_vec;
        proj_ligand_vec_all = input_paras_all.proj_ligand_vec;
        proj_dose_str_vec_all = input_paras_all.proj_dose_str_vec;
        proj_dose_val_vec_all = input_paras_all.proj_dose_val_vec;
        Num_sample_all = input_paras_all.Num_sample;
        save_metric_name = 'Sim12_dual_signle_ligand_codon_metric.mat';
        cal_codon =1;
        ncpu=10;
        pc=parcluster('local');
        pc.NumThreads=2;%
        parpool(pc,ncpu)
        
        parfor i_comb = 1:10
            index_dual_ligand_match = i_comb;%1;2;6;7;11;
            input_paras = struct();
            input_paras.proj_num_vec = proj_num_vec_all(i_comb);
            input_paras.proj_ligand_vec = proj_ligand_vec_all(i_comb);
            input_paras.proj_dose_str_vec = proj_dose_str_vec_all(i_comb);
            input_paras.proj_dose_val_vec = proj_dose_val_vec_all(i_comb);
            input_paras.Num_sample = 1000;%Num_sample_all(i_comb);% 3;%input_paras.Num_sample(index_dual_ligand_match);
            
            [sim_data_tbl{i_comb},data{i_comb},metrics{i_comb},collect_feature_vects{i_comb}] = para_sample_sim_single_2_ligand_2023_12(data_save_file_path,input_paras,cal_codon,'alldose')
        end
        delete(gcp)
        
        save(strcat(data_save_file_path,save_metric_name),'sim_data_tbl','data','metrics','collect_feature_vects','-v7.3');
        
    end
    
    
    % scramble dual ligand stim paramters to be deleted
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % run all conditions
        % Initializing_sampling_supriya_data_2023
        
        Initialize_combinatorial_ligands
        % input_paras.Num_sample = [3];
        save_metric_name = 'Sim11_all_comb_scramble_codon_metric.mat';
        cal_codon =1;
        % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras)
        para_sample_scramble_fitting_sim_sti_doses_var_2023_12(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        
    end
    
    % total scramble dual ligand stim paramters
    
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % run all conditions
        % Initializing_sampling_supriya_data_2023
        
        Initialize_combinatorial_ligands
        % input_paras.Num_sample = [3];
        save_metric_name = 'Sim11_all_comb_all_scramble_codon_metric_202402.mat';
        cal_codon =1;
        % input_paras.Num_sample = 2;
        % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras)
        
        para_sample_all_scramble_fitting_sim_sti_doses_var_2024_02(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    end
    
    % 5 ligand stimu to delete?
    if 0
        monolix_data_save_file_path = '../SAEM_proj_2023/';
        % Initializing_sampling_supriya_data_2023
        input_paras.proj_ligand_vec = {{'TNF' 'LPS' 'CpG' 'PolyIC' 'Pam3CSK'}};
        input_paras.proj_dose_str_vec = {{'10ng/mL','10ng/mL','100nM','100ug/mL','100ng/mL'}};
        input_paras.proj_dose_val_vec = {{10,10,100,100000,100}};
        input_paras.Num_sample = 200;%3 ligand
        input_paras.proj_num_vec = {[2,3,4,5,6]};
        
        save_metric_name = 'Sim7_5ligands_codon_metric.mat';
        cal_codon = 1;
        para_sample_fitting_sim_sti_doses_var_2023_11(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
        
    end
end

%% Figure Dual-ligand pred

%% single ligand sampling for Sjroen syndrome
if 0 % sampling
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
    
    if 0
    
    for i_r =1:2
        
        save_metric_name = strcat('Sim5_SS_codon_metric_p25x_r',num2str(i_r),'.mat');
        para_fitting_sim_IkBao_202401(0.25, data_save_file_path,input_paras,cal_codon,save_metric_name);
        
    end
    end
    
    %
    %     i_name = 1;
    %     for i_r =1:2
    %
    %         % params6 = 0.1 * params6_wt;
    %         save_metric_name{i_name} = strcat('Sim5_SS_codon_metric_p1x_r',num2str(i_r),'.mat');
    %         fold_params6(i_name) = 0.1;
    %         i_name = i_name+1;
    %
    %         % params6 = 0.01 * params6_wt;
    %         save_metric_name{i_name} = strcat('Sim5_SS_codon_metric_p01x_r',num2str(i_r),'.mat');
    %         fold_params6(i_name) = 0.01;
    %         i_name = i_name+1;
    %
    %         % params6 = 0 * params6_wt;
    %         save_metric_name{i_name} = strcat('Sim5_SS_codon_metric_0x_r',num2str(i_r),'.mat');
    %         fold_params6(i_name) = 0;
    %         i_name = i_name+1;
    %     end
    %
    %     ncpu=6;
    %     pc=parcluster('local');
    %     pc.NumThreads=2;%
    %     parpool(pc,ncpu)
    %
    %     parfor i_name = 1:length(save_metric_name)
    %         para_fitting_sim_IkBao_202401(fold_params6(i_name), data_save_file_path,input_paras,cal_codon,save_metric_name{i_name});
    %     end
    %     delete(gcp)
    
end

if 0 %transfer the data into machine learning format
    transfer_data_machine_learning_cal_format_SS_sampling
    
end

if 0 %debug
    addpath('./debug/')
    
    check_heatmap_SS_p1x
end

if 0 % visualize machine learning results
    visualize_classification_ML_SS_sampling
end

%% Application: MI within the network

if 0
    transfer_diff_var_codon_mutual_info_format_0331
    % before 03/31/2024
    % transfer_diff_var_codon_mutual_info_format
end

if 0 %no noise
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

if 0 %only TAK noise Sim22
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

if 0 %only core module noise
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

if 0 %only rcp noise
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

if 0
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
    
    var_input.paranames = {'params52n2','params65n2'};%LPS
    var_input.paravals = [1,1]';
    
    var_input.specienames = {};
    var_input.specievals = [];
    
    
    cal_codon =1;
    
    i_r =1;
    save_metric_name = strcat('Sim20_reduce_TAK1ac_heterogeneity_codon_metric',num2str(i_r),'.mat');
    para_fitting_sim_alter_var_202402(var_input, data_save_file_path,input_paras,cal_codon,save_metric_name);
    
    
end


if 0 %% receptor reduce heterogeneity
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

if 0 %% receptor reduce heterogeneity for IkBas/s
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

if 0 % control: wt, IkB-/-
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

if 0 % generate data, reducing NFkB var, using representaive, vs 9 random picked selected cells
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
    
    ncpu=5;
    pc=parcluster('local');
    pc.NumThreads=2;%
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


% before 01/31/2024
% transfer time trajectory into MI format
if 0
    transfer_sample_traj_into_mutual_info_format
end

if 0 % plot the trajs and codon distribution
    transfer_sample_traj_into_mutual_info_format
end

%% publish: not figures, numbers in the manuscript

%% calculate responder ratio: from supriya's data
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