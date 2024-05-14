function [] = draw_para_distribution_sampling(data_save_file_path,fig_save_path)

% data_proj_num_vec=[18,25:28];%TNF,LPS [27];%


%% para_meter_distribution
%
sti_vec = {'TNF'};%,'LPS','CpG','polyIC','Pam3CSK'};

% sti_vec = {'TNF'};% ,'Pam3CSK'};
proj_num_vec =[18];%,25:28];
para_num = [6,7,6,7,6];

%% sampling
mu =        [0.676, 0.232,      0.0293,     0.245,      4.53e-6,    0.00423 ]; % mean value
sigma =     [1.83,  2.47,       8.19,       2.06,       1.92,       1.91]; %standard deviation of the ramdom effects
corr_mat = [1,      -0.116      -0.234      0.0497      0.381       0.558
            -0.116  1           0.929       0.852       -0.134      0.334
            -0.234  0.929       1           0.846       -0.154      0.182
            0.0497  0.852       0.846       1           -0.0229     0.217
            0.381   -0.134      -0.154      -0.0229     1           0.36
            0.558   0.334       0.182       0.217       0.36        1];
para_min = [0.05188,0.01,0.01,0.01,8.224e-7,0.00116];

para_max = [5.188,1,1,1,8.224e-5,0.116];

sigma_mat = diag(sigma)*corr_mat*diag(sigma);

mu_log = log((mu - para_min)./(para_max-mu));

N_cells = 1500;

R = mvnrnd(mu_log,sigma_mat,N_cells);

para_val = (exp(R).*(ones(N_cells,1)*para_max) + ones(N_cells,1)*para_min)./(1+exp(R));

%% 
for i_sti = 1:length(sti_vec)
    proj_num = proj_num_vec(i_sti);
    [~,para_name,estimates] = read_draw_individual_parameters(proj_num,data_save_file_path);
    plot_est_distribution(para_val, para_name(1:para_num(i_sti)),estimates)
    
    figure(2)
    saveas(gcf, strcat(fig_save_path,sti_vec{i_sti},'_','sampled_parameter_distrib'),'epsc')
    close
    
    
end

