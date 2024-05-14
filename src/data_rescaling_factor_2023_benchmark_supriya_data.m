function rescale_factor = data_rescaling_factor_2023_benchmark_supriya_data(NFkB_max_range,data_path)

% We assume cells in response to LPS 100ng can reach
% the highest Nucleus NFkB concentration (all the NFkB can enter into nucleus).
% (Pam3CSK might not be able to reach the highest Nuc NFkB conc)

%% define the rescale factor: based on highest NFkB response to LPS


% ligand_index = 4; dose_index = 5;
% sti_info = get_stimulus_info(ligand_index,dose_index,params.exp_info);
    
% data_index_all = (dataTbl.Ligand == sti_info.ligand_name);
data = readmatrix(strcat(data_path,'Bench_Supriya_LPS1.xlsx'));
data2 = readmatrix(strcat(data_path,'Bench_Supriya_LPS2.xlsx'));

rescale_factor1 = NFkB_max_range(2) / prctile(max(data,[],2),90) ;% highest dose of LPS
rescale_factor2 = NFkB_max_range(2) / prctile(max(data2,[],2),90) ;% highest dose of LPS
rescale_factor = (rescale_factor1+rescale_factor2)/2;