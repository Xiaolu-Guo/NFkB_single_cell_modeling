function [] = Sim_dual_ligand_stim()
% co-sti simulation, using the original fitted parameters

%_supriya_data_sampling
run_co_sti = 1;
draw_co_sti_sampling = 0;

%% initializing

data_save_file_path = '../raw_data2023/';%_fay_parameter/';
fig_save_path = '../SubFigures2023/';

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path);
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

if run_co_sti
    % run all conditions
    % Initializing_sampling_supriya_data_2023
    
    Initializing_sampling_supriya_data_2023_05
    % single: 1,2,4,6,7
    % dual: 3,5,9,10,14,15,20,26,27
    % tri 21,22
    
    dual_ligand_length = 10;
    single_ligand_length = 5;


   % input_paras.Num_sample = 1000;
    
    index_sim = [3;5;9;10;14;15;20;26;27;29;1;2;4;6;7];%1;2;6;7;11;
    input_paras.proj_num_vec = input_paras.proj_num_vec(index_sim);
    input_paras.proj_ligand_vec = input_paras.proj_ligand_vec(index_sim);
    input_paras.proj_dose_str_vec = input_paras.proj_dose_str_vec(index_sim);
    input_paras.proj_dose_val_vec = input_paras.proj_dose_val_vec(index_sim);
    input_paras.Num_sample = [1*ones(dual_ligand_length,1);
        2*ones(single_ligand_length,1)]*300;
    input_paras.var_fold_mat = 1;
    input_paras.var_fold_vec = 1;
    
    save_metric_name = 'Sim5_codon_all5dose_metric.mat';
    cal_codon =1;
    % para_sample_fitting_sim_sti_doses_var_2023_05(data_save_file_path,monolix_data_save_file_path,input_paras)
    para_sample_fitting_sim_sti_doses_var_2024_05(data_save_file_path,input_paras,cal_codon,save_metric_name,'alldose')
    
    vers_fig = '2023';
    single_or_dual = 'dual';
    %     para_sample_sim(data_save_file_path,Num_sample)
    %     draw_sampling_traj(save_filename,data_save_file_path, vers_fig, fig_save_path)
end

end
