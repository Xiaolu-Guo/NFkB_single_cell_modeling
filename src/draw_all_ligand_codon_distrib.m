function [] = draw_all_ligand_codon_distrib(data_save_file_path,fig_save_path)
%This is for visualizing Monolix estimation
% all; -66,-99,-100,-101,-99-100-101,-receptor1,-rpt2,-rcpt3

fig_opt.paper_opt.paperpos=[0,0,150,180]*3;
fig_opt.paper_opt.papersize=[150 180]*3;

    % data_proj_num_vec = data_proj_nums{i_ligand};    

    load(strcat(data_save_file_path,'All_ligand_codon_2.mat'))% data,collect_feature_vects,metrics));
       
        %% codon distribution
    % violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
    % codon and sti seperately, and save
    fig_opt.save_file = strcat(fig_save_path,'All_ligand_codon_2_distrib');
    violin_plot_codon(collect_feature_vects,fig_opt) % draw all codon and sti in one

 




