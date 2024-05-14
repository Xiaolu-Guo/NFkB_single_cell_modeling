function [] = draw_para_distri_hist_2(data_save_file_path,fig_save_path)

% data_proj_num_vec=[18,25:28];%TNF,LPS [27];%


%% para_meter_distribution
%
sti_vec = {'TNF','LPS','CpG','polyIC','Pam3CSK'};

% sti_vec = {'TNF'};% ,'Pam3CSK'};
proj_num_vec =[3:7];
para_num = [6,7,6,7,6];

cv_mat = zeros(5,7);
prct_width_mat = zeros(5,7);
for i_sti = 1:length(sti_vec)
    proj_num = proj_num_vec(i_sti);
    [para_sample{i_sti},para_name{i_sti},estimates{i_sti}] = read_draw_individual_parameters_2(proj_num,data_save_file_path);
    
    for i_para = 1:para_num(i_sti)
        cv_mat(i_sti,i_para) = std(para_sample{i_sti}(:,i_para))/mean(para_sample{i_sti}(:,i_para));
        prct_width_mat(i_sti,i_para) = (prctile(para_sample{i_sti}(:,i_para),95)-prctile(para_sample{i_sti}(:,i_para),5))...
            /prctile(para_sample{i_sti}(:,i_para),50);
    end
    
    for i_para = 1:para_num(i_sti)
        if 0
            figure(1)
            [~,edges_i] = histcounts(log10(para_sample{i_sti}(:,i_para)));
            paperpos = [0,0,60,50]*1.8;
            papersize = [60,50]*1.8;
            draw_pos = [10,10,40,30]*1.8;
            %histogram(para_sample{i_sti}(:,i_para),10.^edges_i,'Normalization','pdf');hold on
            %histogram(para_sample{i_sti}(:,i_para),10.^edges_i,'Normalization','pdf');hold on
            histogram(para_sample{i_sti}(:,i_para),10.^edges_i,'Normalization','probability');hold on
            
            set(gcf, 'PaperUnits','points')
            set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
            set(gca,'fontsize',6,'fontweight','b','XScale','log')
            % xlabel(para_name{i_sti}{i_para})
            % ylabel('probability')
            saveas(gcf,strcat(fig_save_path,'parameter_distr_',sti_vec{i_sti},'_',para_name{i_sti}{i_para}),'epsc')
            close
        end
        
    end
    
end

sti_vec = {'TNF','LPS','CpG','polyIC','Pam3CSK'};

c = categorical(["TNF" "LPS" "CpG"...
    "polyIC" "Pam3CSK"]);
c = reordercats(c,{'TNF','LPS','CpG','polyIC','Pam3CSK' });
figure(3)
b=bar(c,cv_mat' ,'BaseValue',0,'EdgeColor',[0 0 0],'LineWidth',1);%cv_mat,prct_width_mat
ylabel({'90% percentile',' width'})
ylabel('CV')

% c = categorical(para_name{i_sti}(1:para_num(i_sti)));%["c66" "c99" "c100"... "c101" "r1" "r2" "r3"]);
% c = reordercats(c,para_name{i_sti}(1:para_num(i_sti)));

for i_sti=1:length(sti_vec)
    if 1
    figure(2)
    
    paperpos = [0,0,80,50]*1.8;
    papersize = [80,50]*1.8;
    draw_pos = [10,10,60,30]*1.8;
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    c = categorical(para_name{i_sti}(1:para_num(i_sti)));%["c66" "c99" "c100"... "c101" "r1" "r2" "r3"]);
    c = reordercats(c,para_name{i_sti}(1:para_num(i_sti)));
    b=bar(c,cv_mat(i_sti,1:para_num(i_sti)) ,'BaseValue',0,'EdgeColor',[0 0 0],'LineWidth',1);
        set(gca,'fontsize',6,'fontweight','b')
ylabel('CV')
    saveas(gcf,strcat(fig_save_path,'parameter_var_',sti_vec{i_sti}),'epsc')
    close
    end
    
    if 1
        figure(2)
    
    paperpos = [0,0,80,50]*1.8;
    papersize = [80,50]*1.8;
    draw_pos = [10,10,60,30]*1.8;
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    c = categorical(para_name{i_sti}(1:para_num(i_sti)));%["c66" "c99" "c100"... "c101" "r1" "r2" "r3"]);
    c = reordercats(c,para_name{i_sti}(1:para_num(i_sti)));
    b=bar(c,prct_width_mat(i_sti,1:para_num(i_sti)) ,'BaseValue',0,'EdgeColor',[0 0 0],'LineWidth',1);
        set(gca,'fontsize',6,'fontweight','b')
ylabel({'90% percentile',' width'})
    saveas(gcf,strcat(fig_save_path,'parameter_percentile_width_',sti_vec{i_sti}),'epsc')
    close
    end
end
a =1 ;
% xlimit_vec = {[0.05,5.2];[0.01,1];[0.01,1];[0.01,1]};
%
% for i_hist = 1:4
%     figure(i_hist+2)
%     set(gcf, 'PaperUnits','points')
%     set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
%     xlim(xlimit_vec{i_hist})
%
%     xlabel(para_name{i_hist})
%     ylabel('pdf')
%     set(gca,'fontsize',6,'fontweight','b','XScale','log')
%     if i_hist==1
%         legend(sti_vec)
%     end
%     saveas(gcf,strcat(fig_save_path,'core_parameter_distr_',para_name{i_hist}),'epsc')
%     close
% end

%figure(2)
%     saveas(gcf, strcat(fig_save_path,sti_vec{i_sti},'_','parameter_distrib'),'epsc')
%     close
%
%     for i_hist =1:4
%         figure(i_hist+2)
%         [~,edges_i] = histcounts(log10(para_sample(:,i_hist)));
%         histogram(para_sample(:,i_hist),10.^edges_i,'Normalization','pdf');hold on
%
%     end
%
%
%     data_save_file_path = '../SAEM_proj_2022/';
%     [para_sample,para_name,estimates] = read_draw_individual_parameters_seperate_dose(proj_num,data_save_file_path);
%     stimuli_info_tbl = get_stimuli_info_tbl();
%     estimates.cell_num = stimuli_info_tbl.num_cells((stimuli_info_tbl.Ligand == sti_vec{i_sti}));
%     estimates.dose_color = {'b','g','r','c','m'};%,'b','g','r'};
%     plot_est_distribution_seperate_dose(para_sample(:,1:para_num(i_sti)), para_name(1:para_num(i_sti)),estimates)
%
%     figure(2)
%     saveas(gcf, strcat(fig_save_path,sti_vec{i_sti},'_','parameter_distrib_seperate_dose'),'epsc')
%     close
%
