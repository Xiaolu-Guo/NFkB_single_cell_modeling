% For Guo et al. Figure 2D & Figure S4D, Machine learning results
% resutls from python codes, for random forest classifier using signaling
% codon as input and ligand infomation as output
% 
% tested 05/09/2024, Matlab 2020a

i_data_set = 1;
data_name{i_data_set} = 'Ade';
data_set_all = {'exp','fitting','sampling'};
method_all = {'random_forest','random_forest','random_forest'};

confusion_mat_all{1} = [172   4   5   4  16
    13 116  12  37  15
    5  31  49  28  27
    2  54  12  56  11
    14  26   8  10 113];% random forest

confusion_mat_all{2} = [197   2   2   0   0
    6 129  14  31  13
    5  27  63  22  23
    1  45  12  66  11
    1  13  16   3 138];% random forest Fitting

confusion_mat_all{3} =  [193   5   1   2   1
    11 113  32  30  11
    4  28 102  37  25
    7  25  26 128  15
    1   7  26   9 151];% random forest Sampling

for i_data_set = 1:length(confusion_mat_all)
    confusion_mat = (confusion_mat_all{i_data_set})';
    data_set = data_set_all{i_data_set};
    method = method_all{i_data_set};
    confusion_score = confusion_mat./(sum(confusion_mat,2)*ones(1,size(confusion_mat,2)));
    
    ligand_all = {'TNF','Pam3CSK','CpG', 'LPS','PolyIC'};
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];

    paper_pos = [0,0,150,120]*1.3;
    paper_size = [150,120]*1.3;
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    h = heatmap( ligand_all,ligand_all,confusion_score,'Colormap',mymap,'CellLabelColor','k');%'none'
    caxis([0,1])
    % h.FontColor = [0,0,0];
    
    %         for i_index = 1:length(index)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.8,0.8,0.8}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8} {red}
    %         end
    %
    %         for i_index = 1:length(index_non_wide)
    %             h.XDisplayLabels{i_index} = ['\color[rgb]{0.4,0.4,0.4}' h.XDisplayLabels{i_index}];%[rgb]{0.8,0.8,0.8}
    %         end
    h.CellLabelColor = [0,0,0];
    h.CellLabelFormat = '%0.2g';
    h.FontColor = [0,0,0];

    set(gca,'FontSize',7,'FontName','Arial')
    saveas(gcf,strcat(fig_save_path,'percision_',data_set,'_',method),'epsc')
    close()
end

% fscore(i_ligand) = confusion_mat(

