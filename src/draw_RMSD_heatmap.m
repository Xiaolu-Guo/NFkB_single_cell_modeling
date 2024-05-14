function [] = draw_RMSD_heatmap(data_save_file_path,fig_save_path) 


i_row = 1;
i_row_vec = 1;

vers = '_allsti_large_change';
title_h = '0.1x and 10x';
% small_fold = 0.1;
% large_fold = 10;
plot_num = 2;
% 
% vers = '_allsti_small_change';
% title_h = '0.9x and 1.1x';
% % small_fold = 0.9;
% % large_fold = 1.1;
% plot_num = 1;

data_set = {strcat('TNF',vers,'.mat');
    strcat('CpG',vers,'.mat');
    strcat('Pam3CSK',vers,'.mat');
    strcat('LPS',vers,'.mat');
    strcat('polyIC',vers,'.mat');
    strcat('coreup',vers,'.mat');
    strcat('core',vers,'.mat')};

for i_data_set =1:length(data_set)
    
    load(strcat(data_save_file_path,data_set{i_data_set}));
    % CodewordList = table2cell(FeatureListTable(:,3));
    % codewords = unique(CodewordList);
    
    parameter_vecs = unique(sim_table_data.parameter_name,'stable');
    
    for i_para = 1:length(parameter_vecs)
        index_para = find(strcmp(sim_table_data.parameter_name,parameter_vecs{i_para}));
        ligand_vecs = unique(sim_table_data.ligand(index_para),'stable');
        i_column = 1;
 
        for i_ligand = 1:length(ligand_vecs)
            
            index_ligand = find(strcmp(sim_table_data.parameter_name,parameter_vecs{i_para}) &...
                strcmp(sim_table_data.ligand,ligand_vecs{i_ligand}));
            dose_vecs = unique(sim_table_data.dose_str(index_ligand),'stable');
            
            for i_dose = 1:length(dose_vecs)
                
                index = strcmp(sim_table_data.parameter_name,parameter_vecs{i_para}) &...
                    strcmp(sim_table_data.dose_str,dose_vecs{i_dose}) &...
                    strcmp(sim_table_data.ligand,ligand_vecs{i_ligand});
                
                index0 =  index &...
                    (cell2mat(sim_table_data.parameter_fold_change) == 1);
                
                index_small = index &...
                    (cell2mat(sim_table_data.parameter_fold_change) <1 );% == small_fold);
                
                index_large = index &...
                    (cell2mat(sim_table_data.parameter_fold_change) >1) ;% == large_fold);
                
                
                
                rmsd_1p1 = sqrt(sum((sum(cell2mat(sim_table_data.trajectory(index_large))) -...
                    sum(cell2mat(sim_table_data.trajectory(index0)))).^2/...
                    (length(sum(cell2mat(sim_table_data.trajectory(index0))))-1)));
                rmsd_0p9 = sqrt(sum((sum(cell2mat(sim_table_data.trajectory(index_small))) -...
                    sum(cell2mat(sim_table_data.trajectory(index0)))).^2/...
                    (length(sum(cell2mat(sim_table_data.trajectory(index0))))-1)));
                
                sens_mat(i_row,i_column) = max(rmsd_1p1,rmsd_0p9);
                
                
                if ~(sum(index_small) * sum(index_large))
                    sens_mat(i_row,i_column) = NaN;
                end
                
                
                i_column = i_column+1;
            end
        end
        para_names.index{i_row} = i_row;
        para_names.name{i_row} = parameter_vecs{i_para};
        i_row = i_row+1;
    end
    
    i_row_vec(end+1)  = i_row;
    
end
pick_sens_parameters

sens_mat_draw = sens_mat(~isnan(sens_mat(:,1)),[1:4,6,8:12,14,16:19]);
para_names.index = (1:sum(~isnan(sens_mat(:,1))))';
para_names.name =para_names.name(~isnan(sens_mat(:,1)));
% subplot(1,2,plot_num)
h = heatmap(sens_mat_draw,'ColorMap',parula);

% ticklabels = get(gca,'YTickLabel');
% % prepend a color for each tick label
% ticklabels_new = cell(size(ticklabels));

para_sens = {'params53';'params54';...
    'params85';'params93';...
    'params68';'params75';...
    'params35';'params40';'params44';...
    'params77';'params79';'params83';...
    'params66';'params99';'params100';'params101'};

YLabels = 1:length(para_names.index);
XLabels = 1:15;
% Convert each number in the array into a string
CustomYLabels = string(YLabels);
CustomXLabels = string(XLabels); 

% Replace all but the fifth elements by spaces
% Set the 'XDisplayLabels' property of the heatmap 
% object 'h' to the custom x-axis tick labels
para_name_table = table(para_names.index(:),para_names.name(:),...
    'VariableNames',fieldnames(para_names));

CustomYLabels(:) = " ";
CustomXLabels(:) = " ";


for j_para_sens = 1:length(para_sens)
index_vec = para_names.index(strcmp(para_name_table.name,para_sens{j_para_sens}));
CustomYLabels(index_vec) = para_sens{j_para_sens};

end

sti_vec = {'TNF','CpG','Pam3CSK','LPS','polyIC'};
indexX_vec = [2,5,8,11,14];

for i_sti = 1:length(sti_vec)
  CustomXLabels(indexX_vec(i_sti)) = sti_vec{i_sti};

end  
    
h.YDisplayLabels = CustomYLabels;
h.XDisplayLabels = CustomXLabels;
title(title_h)

end