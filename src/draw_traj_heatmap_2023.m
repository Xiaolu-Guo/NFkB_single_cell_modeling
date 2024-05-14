function [] =  draw_traj_heatmap_2023(data_save_file_path,fig_save_path) 

%This is for visualizing Monolix estimation

%% filepath


% clear all


    % data_save_file_path = '../SAEM_proj_2022_2/';
data_proj_num_vec = [2:6];
% data_proj_num_vec = 5;

% load('allsti_data_codon_0603.mat')

%%
data = read_monolix_data_to_matrix_2023(data_proj_num_vec,data_save_file_path);

vis_data_field = {'exp_mode_filter_nan','pred_mode_filter_nan'};%,'sample'};
data_label = {'experiments','prediction'};%,'sample'};

% [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter

% savepath='./draft_figures/';



%% trajectories
vis_data_field = {'exp','pred_mode'};
% 
% for i_sti = length(data.(vis_data_field{1}))
%  [~, data.order{i_sti}] =sort(collect_feature_vects.Speed{2*i_sti-1},'descend');%OscVsNonOsc PeakAmplitude
% 
% end
% % 
% plot_data(data,vis_data_field)

plot_data_seperate_ligand(data,vis_data_field,fig_save_path)
% cell_num_vec = [1:5;6:10;11:15;16:20;21:25;26:30]
% plot_trajectory(data,vis_data_field,cell_num_vec)
