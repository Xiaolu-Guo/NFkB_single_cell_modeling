function [] = draw_para_distribution(data_save_file_path,fig_save_path)

% data_proj_num_vec=[18,25:28];%TNF,LPS [27];%


%% para_meter_distribution
%
sti_vec = {'TNF','LPS','CpG','polyIC','Pam3CSK'};

% sti_vec = {'TNF'};% ,'Pam3CSK'};
proj_num_vec =[18,25:28];
para_num = [6,7,6,7,6];

for i_sti = 1:length(sti_vec)
    proj_num = proj_num_vec(i_sti);
    [para_sample,para_name,estimates] = read_draw_individual_parameters(proj_num,data_save_file_path);
    plot_est_distribution(para_sample(:,1:para_num(i_sti)), para_name(1:para_num(i_sti)),estimates)
    
    figure(2)
    saveas(gcf, strcat(fig_save_path,sti_vec{i_sti},'_','parameter_distrib'),'epsc')
    close
    
    for i_hist =1:4
        figure(i_hist+2)
        [~,edges_i] = histcounts(log10(para_sample(:,i_hist)));
        histogram(para_sample(:,i_hist),10.^edges_i,'Normalization','pdf');hold on
        
    end
    
    
    data_save_file_path = '../SAEM_proj_2022/';
    [para_sample,para_name,estimates] = read_draw_individual_parameters_seperate_dose(proj_num,data_save_file_path);
    stimuli_info_tbl = get_stimuli_info_tbl();
    estimates.cell_num = stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == sti_vec{i_sti}));
    estimates.dose_color = {'b','g','r','c','m'};%,'b','g','r'};
    plot_est_distribution_seperate_dose(para_sample(:,1:para_num(i_sti)), para_name(1:para_num(i_sti)),estimates)
    
    figure(2)
    saveas(gcf, strcat(fig_save_path,sti_vec{i_sti},'_','parameter_distrib_seperate_dose'),'epsc')
    close
    
end

paperpos = [0,0,60,50]*1.8;
papersize = [60,50]*1.8;
draw_pos = [10,10,40,30]*1.8;
xlimit_vec = {[0.05,5.2];[0.01,1];[0.01,1];[0.01,1]};

for i_hist = 1:4
    figure(i_hist+2)
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    xlim(xlimit_vec{i_hist})
    
    xlabel(para_name{i_hist})
    ylabel('pdf')
    set(gca,'fontsize',6,'fontweight','b','XScale','log')
    if i_hist==1
        legend(sti_vec)
    end
    saveas(gcf,strcat(fig_save_path,'core_parameter_distr_',para_name{i_hist}),'epsc')
    close
end
