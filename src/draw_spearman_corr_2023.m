function [] = draw_spearman_corr_2023(data_save_file_path,fig_save_path)

load(strcat(data_save_file_path,'All_ligand_codon_2023.mat'))% data,collect_feature_vects,metrics));
mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
    ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];


collect_feature_vects.info_ligand_ctg = categorical(collect_feature_vects.info_ligand);
collect_feature_vects.info_dose_str_ctg = categorical(cellfun(@char,collect_feature_vects.info_dose_str,'UniformOutput',false));
collect_feature_vects.info_data_type_ctg = categorical(collect_feature_vects.info_data_type);

data.info_ligand_ctg = categorical(data.info_ligand);
data.info_dose_str_ctg = categorical(cellfun(@char,data.info_dose_str,'UniformOutput',false));

features = {'Duration','EarlyVsLate','OscVsNonOsc','PeakAmplitude','Speed','TotalActivity'};

sti_vec = {'TNF','LPS','CpG','PolyIC','Pam3CSK'};

for i_sti = 1:length(sti_vec)
    ligand = sti_vec{i_sti};
    dose_vec = unique(collect_feature_vects.info_dose_str_ctg(collect_feature_vects.info_ligand_ctg==ligand),'stable');
    for i_dose = 1:length(dose_vec)
        dose  = dose_vec(i_dose);
        para_vals = data.parameters_mode_nan{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose};
        para_vals = para_vals(:,1:end-1);
        
        %% core para
        corr_mat =  corr(para_vals(:,[1:3,size(para_vals,2)]),...
            [collect_feature_vects.(features{1}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{2}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{3}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{4}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{5}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{6}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'}],...
            'Type','Spearman');
        para_names_str = data.para_name{data.info_ligand_ctg == ligand & data.info_dose_str_ctg == dose}(:,1:end-1);
        para_names_str{end} = 'NFKBInit';
        figure(1)
        set(gcf, 'PaperUnits','points')
        set(gcf, 'PaperPosition', [0,0,300,50+4*35],'PaperSize', [300,50+4*35])%,'Position',draw_pos [20,20,280,280]
        h = heatmap( features,para_names_str([1:3,size(para_vals,2)]),corr_mat,'Colormap',mymap,'CellLabelColor','k');%'none'
        h.Title = strcat(ligand,'-',char(dose));
        caxis([-1,1])
        
        saveas(gcf,strcat(fig_save_path,'corr_core_',ligand,'_',replace(char(dose),'/','p')),'epsc')
        close()
        
        %% receptor para
        corr_mat =  corr(para_vals(:,setdiff(1:size(para_vals,2),[1:3,size(para_vals,2)],'stable')),...
            [collect_feature_vects.(features{1}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{2}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{3}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{4}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{5}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'},...
            collect_feature_vects.(features{6}){collect_feature_vects.info_ligand_ctg==ligand & collect_feature_vects.info_dose_str_ctg == dose & collect_feature_vects.info_data_type_ctg == 'prediction'}],...
            'Type','Spearman');
        figure(2)
        set(gcf, 'PaperUnits','points')
        
        paper_pos = [0,0,300,60+(size(para_vals,2)-4)*35];
        paper_size = [300,60+(size(para_vals,2)-4)*35];
        set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
        h = heatmap( features,para_names_str(setdiff(1:size(para_vals,2),[1:3,size(para_vals,2)],'stable')),corr_mat,'Colormap',mymap,'CellLabelColor','k');%'none'
        h.Title = strcat(ligand,'-',char(dose),'-rcp');
        caxis([-1,1])
        
        saveas(gcf,strcat(fig_save_path,'corr_rcp_',ligand,'_',replace(char(dose),'/','p')),'epsc')
        close()
    end
end
