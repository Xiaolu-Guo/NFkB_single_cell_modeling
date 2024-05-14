% clear all

xls_file_save = 0;
unstim_file_save = 0;
SAEM_save = 0;

addpath('../NFkB_common/')
addpath('./lib/')

rescale_data = 1;
rescale_save = 1;
data_path = '/Users/admin/Documents/my document/Postdoc projects/Projects/NFkB_para_estm_project/Experiments/Benchmark_Supriya/';

%% rescaling the data and save
if rescale_data
    % the macrophage cell volume is 6 pl, measured by Stefanie Luecke
    NFkB_max_range = [0.25, 1.75]/6;    
    % We assume cells in response to LPS 100ng can reach
    % the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
    % (Pam3CSK might not be able to reach the highest Nuc NFkB conc) 
    rescale_factor = data_rescaling_factor_2023_benchmark_supriya_data(NFkB_max_range,data_path);    
end

%% non-stim data
% data_rescale(rescale_factor,ligand_index,dose_index,params)
% updated on 11/24/2021
% rescale the data using the rescaling factor
data_file_list = {'Bench_Supriya_CpG1.xlsx'
    'Bench_Supriya_CpG2.xlsx'
    'Bench_Supriya_LPS1.xlsx'
    'Bench_Supriya_LPS2.xlsx'
    'Bench_Supriya_P3K1.xlsx'
    'Bench_Supriya_P3K2.xlsx'
    'Bench_Supriya_PIC1.xlsx'
    'Bench_Supriya_PIC2.xlsx'
    'Bench_Supriya_TNF1.xlsx'
    'Bench_Supriya_TNF2.xlsx'};
write_data_file_list = {'rescaled_Bench_Supriya_CpG1.xlsx'
    'rescaled_Bench_Supriya_CpG2.xlsx'
    'rescaled_Bench_Supriya_LPS1.xlsx'
    'rescaled_Bench_Supriya_LPS2.xlsx'
    'rescaled_Bench_Supriya_P3K1.xlsx'
    'rescaled_Bench_Supriya_P3K2.xlsx'
    'rescaled_Bench_Supriya_PIC1.xlsx'
    'rescaled_Bench_Supriya_PIC2.xlsx'
    'rescaled_Bench_Supriya_TNF1.xlsx'
    'rescaled_Bench_Supriya_TNF2.xlsx'};

if rescale_save    
    for i_data_file = 1:length(data_file_list)
        data = readmatrix(strcat(data_path,data_file_list{i_data_file}));       
        % for debugging
        %histogram(max(data,[],2));hold on        
        data_rescale = data * rescale_factor;
        writematrix(data_rescale,strcat(data_path,write_data_file_list{i_data_file}))
        % for debugging
        % histogram(max(data_rescale,[],2));hold on        
    end
end
