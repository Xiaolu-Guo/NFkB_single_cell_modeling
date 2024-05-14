% clear all
% load the codon for all sampled data and supriya and ade's data for
% comparison
% comapring ade's data and supriya's data: the batch effects.
% comparing surpiya's data and model prediction 
% check each chunk to set 1 or 0 to run or not run the corresponding part
% 
cal_codon = 1;

addpath('./lib/')
addpath('./src/')
addpath('./bin/')

fig_save_path = '../SubFigures2023/';

data_save_file_path_1 = '../raw_data2023/';%_fay_parameter/';

%%  calculate/load the codon for all data
if cal_codon
    
    % Ade's data:
    load(strcat(data_save_file_path_1,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
    data_ade = data;
    metrics_ade = metrics;
    
    % Supriya's data
    load('Supriya_data_metrics.mat')
    data_supriya = data;
    metrics_supriya = metrics;
    
    % TNF 10ng/mL       supriya: 03 index: 4;       ade index: 3    sample: 13
    % Pam3CSK 100ng/mL  supriya: 24 index: 27;      ade index: 18   sample: 12
    % CpG 100nM         supriya: 23 index: 22;      ade index: 11   sample: 11
    % LPS 10ng/mL       supriya: 04 index: 6;       ade index: 6    sample: 14
    % PolyIC 100ug/mL   supriya: 20 index: 15;      ade index: 16   sample: 15
    
    % goal: supriya's data, ade's data, fitting data, and sampled data for
    % benchamarking
    % compare the prediction of smapled data with supriya's data
    
    % TNF 10ng/mL Pam3CSK 100ng/mL supriya: 03 index: 3
    % TNF 10ng/mL CpG 100nM supriya: 03 index: 5
    % TNF 10ng/mL LPS 10ng/mL supriya: 04 index: 8
    % TNF 10ng/mL PolyIC 100ug/mL supriya: 04 index: 9
    % Pam3CSK 100ng/mL CpG 100nM supriya: does not apply
    % Pam3CSK 100ng/mL LPS 10ng/mL supriya: 23 index: 26
    % Pam3CSK 100ng/mL PolyIC 100ug/mL supriya: 13 index: 14
    % CpG 100nM LPS 10ng/mL supriya: 23 index: 25
    % CpG 100nM PolyIC 100ug/mL supriya: 21 index: 18
    % LPS 10ng/mL PolyIC 100ug/mL supriya: 21 index: 19
    
    ade_index = [3;18;11;6;16];
    supriya_index = [4;27;22;6;15];
    supriya_dual_ligand_index = [3;5;8;9;26;14;25;18;19];
    
    % sampled data
    data_filename = 'Sim5_codon_all5dose_metric.mat';
    load(strcat(data_save_file_path_1,data_filename))
    sample_dual_ligand_index = [1;2;3;4;10;9;6;8;5;7];
    sample_index = [13;12;11;14;15];
    
    metric_fields_name = fieldnames(metrics{1});
    for i_metric = 1:length(metrics)
        
        for i_metric_fields = 1:length(metric_fields_name)
            metrics_sample{i_metric}.(metric_fields_name{i_metric_fields}) = metrics{i_metric}.(metric_fields_name{i_metric_fields})(1:9:end,:);
        end
    end
    
    data_sample = data;
    % metrics_sample = metrics;
    
    metrics = [metrics_supriya(supriya_index),... Supriya's single liagand data
        metrics_ade(ade_index*2-1),... Ade's single liagand data
        metrics_ade(ade_index*2),... Fitting Ade single liagand data
        metrics_sample(sample_index),... Sampling single liagand data
        metrics_supriya(supriya_dual_ligand_index),... Supriya's dual ligand data
        metrics_sample(sample_dual_ligand_index)]; % sampling dual ligand data
    
    
    data_info.info_ligand = [data_supriya.info_ligand(supriya_index),... Supriya's single liagand data
        reshape(data_ade.info_ligand(ade_index),1,[]),... Ade's single liagand data
        reshape(data_ade.info_ligand(ade_index),1,[]),... Fitting Ade single liagand data
        data_sample.info_ligand(sample_index),... Sampling single liagand data
        data_supriya.info_ligand(supriya_dual_ligand_index),... Supriya's dual ligand data
        data_sample.info_ligand(sample_dual_ligand_index)]; % sampling dual ligand data
    
    
    for i_info_ligand = 1:length(data_info.info_ligand)
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'polyIC','PolyIC');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'Pam3CSK_TNF','TNF_Pam3CSK');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_TNF','TNF_CpG');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_TNF','TNF_LPS');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_TNF','TNF_PolyIC');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_Pam3CSK','Pam3CSK_CpG');  
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_Pam3CSK','Pam3CSK_LPS');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_Pam3CSK','Pam3CSK_PolyIC');
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'LPS_CpG','CpG_LPS');      
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'PolyIC_CpG','CpG_PolyIC');      
        data_info.info_ligand{i_info_ligand} = replace(data_info.info_ligand{i_info_ligand},'CpG_PolyIC','LPS_PolyIC');
     end
    
    
    data_info.info_dose_str = [data_supriya.info_dose_str(supriya_index),... Supriya's single liagand data
        reshape(data_ade.info_dose_str(ade_index),1,[]),... Ade's single liagand data
        reshape(data_ade.info_dose_str(ade_index),1,[]),... Fitting Ade single liagand data
        data_sample.info_dose_str(sample_index),... Sampling single liagand data
        data_supriya.info_dose_str(supriya_dual_ligand_index),... Supriya's dual ligand data
        data_sample.info_dose_str(sample_dual_ligand_index)]; % sampling dual ligand data
    
    data_info.data_label(1:length(supriya_index)) = {'Exp_Supriya'};
    data_info.data_label(end+1:end+length(ade_index)) = {'Exp_Ade'};
    data_info.data_label(end+1:end+length(ade_index)) = {'Fitting'};
    data_info.data_label(end+1:end+length(sample_index)) = {'Sampling'};
    data_info.data_label(end+1:end+length(supriya_dual_ligand_index)) = {'Exp_Supriya_dual_ligand'};
    data_info.data_label(end+1:end+length(sample_dual_ligand_index)) = {'Sampling_dual_ligand'};
    
    
    [collect_feature_vects,metrics_new] = calculate_codon_from_metric2023(data_info,metrics); %,  parameter
    metrics = metrics_new;
    %
    collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
    collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
    collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);
    %
    save(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'),'collect_feature_vects','metrics','data_info')
    % goal: visualize supriya, ade, fitting, sampling (2 seprate plots) as
    % benchmarking; then visualize supriya dual ligand vs sampling.
else
    load(strcat(data_save_file_path_1,'All_codon_dual_ligand.mat'))
end

%% different data set index setting
index_supriya = 1:5;
index_ade = 6:10;
index_fitting = 11:15;
index_sampling = 16:20;
index_supriya_dual = 21:29;
index_sampling_dual = [30:33,35:39];

%% violin plot of single ligand supriya ade and fitting
if 0
    clear index_violin collect_feature_vects_draw
    codon_fields = fieldnames(collect_feature_vects);
    index_violin(1:3:13) = index_supriya; % supriya
    index_violin(2:3:14) = index_ade; % ade
    index_violin(3:3:15) = index_fitting; % fitting
    for i_codon_fields = 1:length(codon_fields)
        collect_feature_vects_draw.(codon_fields{i_codon_fields}) = collect_feature_vects.(codon_fields{i_codon_fields})(index_violin);
    end
    
    fig_opt.save_file = strcat(fig_save_path,'Single_ligand_exp_supriya_ade_fitting_2023_05_distrib');
    fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
    fig_opt.paper_opt.papersize=[220 180]*3;
    
    % violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
    
    violin_plot_codon_sampling_2023_04(collect_feature_vects_draw,fig_opt) % draw all codon and sti in one
    
    
end

%% violin plot of single ligand supriya ade and sampling 
if 0
    clear index_violin collect_feature_vects_draw
    codon_fields = fieldnames(collect_feature_vects);
    index_violin(1:3:13) = index_supriya; % supriya
    index_violin(2:3:14) = index_ade; % ade
    index_violin(3:3:15) = index_sampling; % fitting
    for i_codon_fields = 1:length(codon_fields)
        collect_feature_vects_draw.(codon_fields{i_codon_fields}) = collect_feature_vects.(codon_fields{i_codon_fields})(index_violin);
    end
    
    fig_opt.save_file = strcat(fig_save_path,'Single_ligand_exp_supriya_ade_sampling_2023_05_distrib');
    fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
    fig_opt.paper_opt.papersize=[220 180]*3;
    
    % violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
    
    violin_plot_codon_sampling_2023_04(collect_feature_vects_draw,fig_opt) % draw all codon and sti in one
end

%% violin plot of dual ligands supriya and sampling
if 0
    clear index_violin collect_feature_vects_draw
    codon_fields = fieldnames(collect_feature_vects);
    index_violin(1:2:17) = index_supriya_dual; % supriya dual ligand
    index_violin(2:2:18) = index_sampling_dual; % sampling dual ligand
    
    for i_codon_fields = 1:length(codon_fields)
        collect_feature_vects_draw.(codon_fields{i_codon_fields}) = collect_feature_vects.(codon_fields{i_codon_fields})(index_violin);
    end
    
    fig_opt.save_file = strcat(fig_save_path,'Dual_ligand_exp_supriya_sampling_2023_05_distrib');
    fig_opt.paper_opt.paperpos=[0,0,220,180]*3;
    fig_opt.paper_opt.papersize=[220 180]*3;
    
    % violin_plot_codon_2023(collect_feature_vects,fig_opt) % draw all codon and sti in one
    
    violin_plot_codon_dual_ligand_sampling_2023_05(collect_feature_vects_draw,fig_opt) % draw all codon and sti in one
end

%% codon comparison between supriya's data and ade's data
if 0
    comparing_case = 'Single_ligand_exp_supriya_ade_2023_05_';
    index_data1 = index_supriya;
    index_data2 = index_ade;
    
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    diff_method = {'L2','L1','JSD'};
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    
    TNF_color = [119 180 202]/255;
    LPS_color = [222 78 66]/255;
    CpG_color = [137 180 66]/255;
    PolyIC_cclor = [101 77 123]/255;
    P3CSK_color = [229 129 56]/255;
    
    
    for i_method =1:3
        
        [good_fit_pop,diff_distr,sti_vec] = cal_codon_distri_diff(collect_feature_vects,index_data1,index_data2,diff_method{i_method});
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        paper_pos = [0,0,450,200];
        paper_size = [450,200];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,diff_distr,'Colormap',mymap,'CellLabelColor','none');%'none'
        % caxis([0,0.5])
        
        
        saveas(gcf,strcat(fig_save_path,comparing_case,'_distrib_codon_',diff_method{i_method},'_diff'),'epsc')
        close()
        
        figure(1)
        paperpos = [0,0,120,50]*1.8;
        papersize = [120,50]*1.8;
        draw_pos = [10,10,100,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        c = categorical(collect_feature_vects.info_ligand(index_data1));
        c = reordercats(c,collect_feature_vects.info_ligand(index_data1));
        bar_pct = good_fit_pop;
        b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        ylabel({'discrepancy between' 'supriya''s and ade''s data'})
        XL = xlim();
        % ylim([0,0.5]);
        % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)
        b(1).FaceColor = 'flat';
        
        b(1).CData = [TNF_color;LPS_color;CpG_color;PolyIC_cclor;P3CSK_color];
        
        set(gca,'fontsize',6,'fontweight','b')
        saveas(gcf,strcat(fig_save_path,comparing_case,'_bar_codon_',diff_method{i_method},'_good_fit_pop'),'epsc')
        close
    end
end

%% codon comparison between %% supriya's data and sampling
if 1
    comparing_case = 'Single_ligand_exp_supriya_sampling_2023_05_';
    index_data1 = index_supriya;
    index_data2 = index_sampling;
    
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    diff_method = {'L2','L1','JSD'};
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    
    TNF_color = [119 180 202]/255;
    LPS_color = [222 78 66]/255;
    CpG_color = [137 180 66]/255;
    PolyIC_cclor = [101 77 123]/255;
    P3CSK_color = [229 129 56]/255;
    
    
    for i_method =1:3
        
        [good_fit_pop,diff_distr,sti_vec] = cal_codon_distri_diff(collect_feature_vects,index_data1,index_data2,diff_method{i_method});
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        paper_pos = [0,0,450,200];
        paper_size = [450,200];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,diff_distr,'Colormap',mymap,'CellLabelColor','none');%'none'
        % caxis([0,0.5])
        
        
        saveas(gcf,strcat(fig_save_path,comparing_case,'_distrib_codon_',diff_method{i_method},'_diff'),'epsc')
        close()
        
        figure(1)
        paperpos = [0,0,120,50]*1.8;
        papersize = [120,50]*1.8;
        draw_pos = [10,10,100,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        c = categorical(collect_feature_vects.info_ligand(index_data1));
        c = reordercats(c,collect_feature_vects.info_ligand(index_data1));
        bar_pct = good_fit_pop;
        b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        ylabel({'discrepancy between' 'supriya''s and sampling data'})
        XL = xlim();
        % ylim([0,0.5]);
        % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)
        b(1).FaceColor = 'flat';
        
        b(1).CData = [TNF_color;LPS_color;CpG_color;PolyIC_cclor;P3CSK_color];
        
        set(gca,'fontsize',6,'fontweight','b')
        saveas(gcf,strcat(fig_save_path,comparing_case,'_bar_codon_',diff_method{i_method},'_good_fit_pop'),'epsc')
        close
    end
end

%% codon comparison between  dual ligand supriya's data and sampling
if 1 %% 
    comparing_case = 'Dual_ligand_exp_supriya_sampling_2023_05_';
    index_data1 = index_supriya_dual;
    index_data2 = index_sampling_dual;
    
    mymap = [(0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    diff_method = {'L2','L1','JSD'};
    codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
    
    TNF_color = [119 180 202]/255;
    LPS_color = [222 78 66]/255;
    CpG_color = [137 180 66]/255;
    PolyIC_cclor = [101 77 123]/255;
    P3CSK_color = [229 129 56]/255;
     
    for i_method =1:3
        
        [good_fit_pop,diff_distr,sti_vec] = cal_codon_distri_diff(collect_feature_vects,index_data1,index_data2,diff_method{i_method});
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        paper_pos = [0,0,450,200];
        paper_size = [450,200];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec,codon_list,diff_distr,'Colormap',mymap,'CellLabelColor','none');%'none'
        % caxis([0,0.5])

        saveas(gcf,strcat(fig_save_path,comparing_case,'_distrib_codon_',diff_method{i_method},'_diff'),'epsc')
        close()
        
        figure(1)
        paperpos = [0,0,120,50]*1.8;
        papersize = [120,50]*1.8;
        draw_pos = [10,10,100,30]*1.8;
        
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
        c = categorical(collect_feature_vects.info_ligand(index_data1));
        c = reordercats(c,collect_feature_vects.info_ligand(index_data1));
        bar_pct = good_fit_pop;
        b = bar(c, bar_pct, 'EdgeColor',[0 0 0],'LineWidth',0.5);hold on
        ax2 = gca;
        % ytickformat(ax2, '%g%%');
        % ylim([75,100])
        ylabel({'discrepancy between' 'exp and sampling dual-ligand'})
        XL = xlim();
        % ylim([0,0.5]);
        % plot([b(1).XData(1)-0.5,b(3).XData(5)+0.5],[0.1,0.1],'--r','LineWidth',1.5)
        b(1).FaceColor = 'flat';
        
        % b(1).CData = [TNF_color;LPS_color;CpG_color;PolyIC_cclor;P3CSK_color];
        
        set(gca,'fontsize',6,'fontweight','b')
        saveas(gcf,strcat(fig_save_path,comparing_case,'_bar_codon_',diff_method{i_method},'_good_fit_pop'),'epsc')
        close
    end
end


