
function [] = visualize_para_codon_feature_importance_2023_06(data_save_file_path,data_save_path,fig_save_path)
parper_size_codon_para = [100,100];
parper_size_accuracy = [150,100];
load(strcat(data_save_file_path,'All_ligand_codon_2023_t33_cv_filtered_TNF.mat'))

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};
dose_vec = {{'100pg/mL';'1ng/mL';'10ng/mL'};
    {'10ng/mL';'100ng/mL';'1ug/mL'};
    {'33nM';'100nM';'333nM'};%'10nM';;'1uM'
    {'1ng/mL';'3ng/mL';'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'10ug/mL';'33ug/mL';'100ug/mL'}};
dose_symbol = {'Low','Med','High'};

dose_vec = {{'10ng/mL'};
    {'1ug/mL'};
    {'333nM'};%'10nM';;'1uM'
    {'10ng/mL'};%;'33ng/mL';'100ng/mL'
    {'100ug/mL'}};
dose_symbol = {'High'};

data_dose_str = cellfun(@char,data.info_dose_str,'UniformOutput',false);

collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

mymap = [ ones(21,1),(1:-0.05:0)',(1:-0.05:0)'];


features = {'Speed','PeakAmplitude','Duration','TotalActivity','EarlyVsLate','OscVsNonOsc'};

ligand_vec = {'TNF','Pam3CSK','CpG','LPS','PolyIC'};

stat_methods_vec = {'LinearRegression','RamdomForestRegression','XGBRegression'};

ColorScaling_ = {'log','scaled','scaled'};
c_axis_limit = {log([1e-3,1e5]),[0,0.8],[0.1,0.8]};

para_names = {'rcp syn','C deg','end','TAK ac', 'time d1','time d2','Tot NFkB'};

for i_stat_methods = 1:length(stat_methods_vec)
    stat_methods = stat_methods_vec{i_stat_methods};
    score_ = NaN(length(features),5*length(dose_symbol));
    rand_feature_imp = NaN(length(features),5*length(dose_symbol));
    for i_feature = 1:length(features)
        corr_mat_codon = NaN(7,5*length(dose_symbol));
        i_cond = 1;
        for i_sti = 1:length(ligand_vec)
            ligand = ligand_vec{i_sti};
            % dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
            for i_dose = 1:length(dose_vec{i_sti})
                dose  = dose_vec{i_sti}{i_dose};
                dose_sym = dose_symbol{i_dose};
                
                A = readmatrix(strcat(data_save_path,stat_methods,'_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
                
                switch stat_methods
                    case 'XGBRegression'
                        A_rand = readmatrix(strcat(data_save_path,stat_methods,'_rand_',features{i_feature},'_',ligand,'_',dose_sym,'.csv'));
                        rand_feature_imp(i_feature,i_cond) = mean(A_rand(length(A)+1:end));
                end
                score_(i_feature,i_cond) = A(1);
                len_rcp = length(A)-5;
                corr_mat_codon(4:end,i_cond) = A(2:5);%core
                corr_mat_codon(1:len_rcp,i_cond) = A(6:5+len_rcp);%receptor
                switch stat_methods
                    case 'LinearRegression'
                        corr_mat_codon = abs(corr_mat_codon);
                end
                sti_vec{i_cond} = strcat(ligand,'-',dose_sym);
                sti_vec2{i_cond} = strcat(ligand);
                
                i_cond = i_cond+1;
            end
        end
        
        
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        %         paper_pos = [0,0,450,200];
        %         paper_size = [450,200];
        paper_pos = [0,0,parper_size_codon_para];
        paper_size = parper_size_codon_para;
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        
        h = heatmap( sti_vec2,para_names,corr_mat_codon,'Colormap',mymap,'ColorLimits',c_axis_limit{i_stat_methods},...
            'ColorScaling',ColorScaling_{i_stat_methods},...
            'CellLabelColor','none');%'k' % corr_mat_codon
        
        switch features{i_feature}
            case 'OscVsNonOsc'
                for i_index = 1:length(h.XDisplayLabels)
                    h.XDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
            otherwise
                for i_index = 1:length(h.XDisplayLabels)
                    h.XDisplayLabels{i_index} = ['\color[rgb]{1,1,1}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
        end
        for i_index = 1:length(h.YDisplayLabels)
            h.YDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.YDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
        end
        % xlabel(features{i_feature})
        %         if i_stat_methods>1
        %             caxis(c_axis_limit{i_stat_methods})
        %         end
        %h.CellLabelColor = [0,0,0];
        h.CellLabelFormat = '%0.2g';
        set(gca,'fontsize',7,'fontname','Arial');
        colorbar('off')
        
        saveas(gcf,strcat(fig_save_path,stat_methods,'_',features{i_feature},'_para_feature_importance'),'epsc')
        close()
    end
    
    if 1
        figure(1)
        set(gcf, 'PaperUnits','points')
        
        paper_pos = [0,0,parper_size_accuracy];
        paper_size = parper_size_accuracy;
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        h = heatmap( sti_vec2,features,score_,'Colormap',mymap,'CellLabelColor','k');%'k' % corr_mat_codon
        set(gca,'fontsize',7,'fontname','Arial');
        
        title('Accuracy Score')
        caxis([0,1])
        %h.CellLabelColor = [0,0,0];
        h.CellLabelFormat = '%0.2g';
        for i_index = 1:length(h.XDisplayLabels)
            h.XDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
        end
        
        switch stat_methods
            case 'LinearRegression'
                for i_index = 1:length(h.YDisplayLabels)
                    h.YDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.YDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
            otherwise
                for i_index = 1:length(h.YDisplayLabels)
                    h.YDisplayLabels{i_index} = ['\color[rgb]{1,1,1}' h.YDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
        end
        saveas(gcf,strcat(fig_save_path,stat_methods,'_Accuracy_Score_',features{i_feature},'_para_feature_importance'),'epsc')
        close()
        
        switch stat_methods
            case 'XGBRegression'
                figure(3)
                set(gcf, 'PaperUnits','points')
                
                paper_pos = [0,0,parper_size_accuracy];
                paper_size = parper_size_accuracy;
                set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
                h = heatmap( sti_vec2,features,rand_feature_imp,'Colormap',mymap,'CellLabelColor','k');%'k' % corr_mat_codon
                set(gca,'fontsize',7,'fontname','Arial');
                
                title('rand variable feature importance')
                caxis(c_axis_limit{i_stat_methods})
                %h.CellLabelColor = [0,0,0];
                h.CellLabelFormat = '%0.2g';
                for i_index = 1:length(h.XDisplayLabels)
                    h.XDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
                for i_index = 1:length(h.YDisplayLabels)
                    h.YDisplayLabels{i_index} = ['\color[rgb]{0,0,0}' h.YDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
                end
                saveas(gcf,strcat(fig_save_path,stat_methods,'_rand_variable_importance_',features{i_feature},'_para_feature_importance'),'epsc')
                close()
        end
    end
    
    
end

