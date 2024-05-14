
function [] = draw_paraDistrib_by_clusterNo_2024(data_save_file_path)
monolix_data_save_file_path = '../SAEM_proj_2023/';

vers = '';
load(strcat('../raw_data2023/Sim15_5_signle_ligand_codon_metric',vers,'.mat'))

for epsilon = 1:2
    clusterResults = epsilon_cluster_data(data,metrics,epsilon);
    [cluster_sort,cluster_sort_ind] = sort(clusterResults);
    
    % cluster_3cat = (clusterResults>3)+(clusterResults>2);
    cluster_low = clusterResults<4;
    cluster_med = clusterResults==4;
    cluster_high = clusterResults>4;
    
    %!!!!
    % need to be fixed: parameter NFkB abundance is not recorded.
    para_cluster = sim_data_tbl.parameter_value(1:9:9000,:);
    para_cluster_lowSRS = para_cluster(cluster_low,:);
    para_cluster_medSRS = para_cluster(cluster_med,:);
    para_cluster_highSRS = para_cluster(cluster_high,:);
    idnex_vec = {[1,2,3],[1:12,14:16]};
    module_name = {'core','all'};
    parameter_name_vec = {{'params52','parmas99','params101'},...
        {'params52','parmas99','params101',...
        'params54','parmas58',...
        'params35','parmas44','params36',...
        'params85','parmas93','params88',...
        'params77','params79',...'parmas83'
        'params68','parmas75'}};
    parameter_bioname = {{'TAK-act','Time-d1','Time-d2'},...
        {'TAK-act','Time-d1','Time-d2',...
        'Rcp-Syn','Cplx-deg',...
        'Rcp-Syn','Cplx-deg','End-trans',...
        'Rcp-Syn','Cplx-deg','End-trans',...
        'Rcp-Syn','End-trans',...'Cplx-deg',
        'Rcp-Syn','Cplx-deg'}};
    
    minval_vec = {[0.1,0.01,0.01],...
        [0.1,0.01,0.01,...
        8.224e-7,0.0125,...%TNF
        5.25e-6,0.0012,0.0065681,...%LPS
        2e-7,0.00016,0.0015,...%CpG
        3e-7,0.004,...%PolyIC 7e-5,
        1e-7,0.0004]};
    maxval_vec = {[10,1,1],...
        [10,1,1,...
        8.224e-5,1.25,...%TNF
        0.000525,0.12,0.65681...%LPS
        2e-5,0.016,0.15,...%CpG
        3e-5,0.4,...%PolyIC %0.007,
        1e-5,0.04]};
    meanval_vec = {[0.873444153573564,0.114187915274479,0.312262009561515],...
        [0.873444153573564,0.114187915274479,0.312262009561515,...
        1.15878192481335e-5,1.23438713904099,...%TNF
        5.39786309680691e-5,0.0322456274251284,0.106546902388348...%LPS
        9.58224188556514e-7,0.00155094376959101,0.13543439243563...%CpG
        3.00069208724329e-7,0.277971665801912,...%PolyIC 0.00699999999999978,
        6.12006752256586e-7,0.0142716976575488]};
    sigmaval_vec = {[1.97492145569698,2.7399354237761,1.6503361183589],...
        [1.97492145569698,2.7399354237761,1.6503361183589,...
        2.84516310162448,8.32296038511482,...%TNF
        1.46828423110844,4.0242288727563,2.11254679282708...%LPS
        1.31555832593251,4.47214539974379,5.07863699256483...%CpG
        6.22561704078366,5.302214857795,...%PolyIC 79.2801447033656,
        0.638998945588272,3.95386398060705]};
    
    
    for i_index = 1 :length(idnex_vec)
        para_index = idnex_vec{i_index};% core module: 52, 99, 101
        minval = minval_vec{i_index};
        maxval = maxval_vec{i_index};
        
        para_data = {minmax_trans(para_cluster_lowSRS(:,para_index),minval,maxval),...
            minmax_trans(para_cluster_medSRS(:,para_index),minval,maxval),...
            minmax_trans(para_cluster_highSRS(:,para_index),minval,maxval)};
        
        cal_flag = 'zscore';
        meanval = meanval_vec{i_index};
        sigmaval = sigmaval_vec{i_index};
        para_data = {logit_trans(para_cluster_lowSRS(:,para_index),minval,maxval,meanval,sigmaval,cal_flag),...
            logit_trans(para_cluster_medSRS(:,para_index),minval,maxval,meanval,sigmaval,cal_flag),...
            logit_trans(para_cluster_highSRS(:,para_index),minval,maxval,meanval,sigmaval,cal_flag)};
        
        
        dataclass = {'lowSRS','medSRS','highSRS'};
        color_vec ={[0,0,1],[0,1,0],[1,0,0]};
        xylabel = parameter_bioname{i_index};
        plot_scatter_hist(para_data,color_vec,xylabel,dataclass,strcat('para_distrib_by_cluster_',module_name{i_index},'_epsilon',num2str(epsilon)))
    end
    
    minval = minval_vec{2};%Pam3CSK
    maxval = maxval_vec{2};%Pam3CSK
    meanval = meanval_vec{2};%Pam3CSK
    sigmaval = sigmaval_vec{2};%Pam3CSK
    cal_flag = 'zscore';
    para_index = idnex_vec{2};% core module: 52, 99, 101
    para_cluster_logit = logit_trans(para_cluster(:,para_index),minval,maxval,meanval,sigmaval,cal_flag);
    finite_index = sum(isinf(para_cluster_logit),2)==0;
    clusterResults_finite = clusterResults(finite_index);
   cluster_low_finite = clusterResults_finite<4;
    cluster_med_finite = clusterResults_finite==4;
    cluster_high_finite = clusterResults_finite>4;

    [~,para_pca] = pca(para_cluster_logit(finite_index,:));
    
    para_cluster_lowSRS = para_pca(cluster_low_finite,:);
    para_cluster_medSRS = para_pca(cluster_med_finite,:);
    para_cluster_highSRS = para_pca(cluster_high_finite,:);

    
    figure(2)
    scatter(para_cluster_lowSRS(:,1),para_cluster_lowSRS(:,2),2,'b','filled');hold on
    scatter(para_cluster_medSRS(:,1),para_cluster_medSRS(:,2),2,'g','filled');hold on
    scatter(para_cluster_highSRS(:,1),para_cluster_highSRS(:,2),2,'r','filled');
    xlabel('PC1');
    ylabel('PC2');
    paperpos = [0,0,100,100]*2;
    papersize = [100,100]*2;
    draw_pos = [20,20,60,60]*2;
    set(gcf, 'PaperUnits','points')
    set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
    fig_save_path = '../SubFigures2023/';
    
    set(gca, 'fontsize',7,'XTick',{},'YTick',{})
    saveas(gcf,strcat(fig_save_path,'para_distrib_by_cluster_PCA_epsilon',num2str(epsilon)),'epsc')
    close()
    
end
end

function minmax_para = minmax_trans(para_mat,minval,maxval)
for i_para = 1:size(para_mat,2)
    minmax_para(:,i_para) = ((para_mat(:,i_para)-minval(i_para))./(maxval(i_para) - para_mat(:,i_para)));
end
end

function [logit_para] = logit_trans(para_mat,minval,maxval,meanval,sigmaval,cal_flag)
for i_para = 1:size(para_mat,2)
    logit_para(:,i_para) = log((para_mat(:,i_para)-minval(i_para))./(maxval(i_para) - para_mat(:,i_para)));
end

switch cal_flag
    case 'zscore'
        for i_para = 1:size(para_mat,2)
            logit_para(:,i_para) = (logit_para(:,i_para)...
                - log(meanval(i_para)-minval(i_para))./(maxval(i_para) -meanval(i_para)))...
                /(4*sigmaval(i_para))+0.5;
        end
    case other
        
end
end

function plot_scatter_hist(data,color_vec,xylabel,dataclass,fig_save_name)
% data{1},data{2},data{3},.... must has same columns
% corresponding to color_vec{1}, color_vec{2}, ...
% label as xylabel{1},{2},....
% legend as dataclass{1},dataclass{2},...

% dataclass = {'Low','Medium','High'};
% xylabel = {'TotalActivity','Duration','EarlyVsLate','Speed','PeakAmplitude','OscVsNonOsc'};
% xylabel = {'params1','params2'};

% fig_save_name = 'para_distrib_by_cluster';

h_vec = [];
for i_data =  1:length(data)%[2,3]%
    num_cells = length(data{i_data});
    x_val = -ones(length(xylabel)^2*num_cells,1);
    y_val = x_val;
    for i_xylabel =1:length(xylabel)
        for j_xylabel = 1:length(xylabel)
            x_start = (i_xylabel-1)*6*num_cells + (j_xylabel-1) * num_cells + 1;
            y_start = (j_xylabel-1)*6*num_cells + (i_xylabel-1) * num_cells + 1;
            
            if i_xylabel == j_xylabel
                x_val(x_start:x_start+num_cells - 1 ) = 0;
                y_val(y_start:y_start+num_cells - 1)= 0;
            else
                x_val(x_start:x_start+num_cells - 1) = data{i_data}(:,i_xylabel)+(i_xylabel-1)*1.2;
                y_val(y_start:y_start+num_cells - 1) = data{i_data}(:,i_xylabel)+(i_xylabel-1)*1.2;
            end
        end
        
        codon_data = data{i_data}(:,i_xylabel)+(i_xylabel-1)*1.2;
        pts = linspace(min(min(codon_data),(i_xylabel-1)*1.2),max(max(codon_data),i_xylabel*1.2),100);
        [f_para(i_data,:,i_xylabel),xi_para(i_data,:,i_xylabel)] = ksdensity(codon_data,pts,...
            'Function','pdf','Bandwidth',0.05);%,'Bandwidth',bw
        
    end
    h(i_data) = scatter(x_val,y_val,1,color_vec{i_data},'filled');hold on
    h_vec(end+1) = h(i_data);
end

scale_factor = max(f_para,[],[1 2]);

for i_data = 1:length(data)
    
    for i_xylabel = 1:length(xylabel)
        plot(xi_para(i_data,:,i_xylabel),f_para(i_data,:,i_xylabel)/scale_factor(i_xylabel)+(i_xylabel-1)*1.2,'Color',color_vec{i_data},'LineWidth',1); hold on
    end
end

for i_x = 0:length(xylabel)
    plot([i_x*1.2-0.1,i_x*1.2-0.1],[-0.2,length(xylabel)*1.2],'Color',[0.7,0.7,0.7]);hold on
    plot([-0.2,length(xylabel)*1.2],[i_x*1.2-0.1,i_x*1.2-0.1],'Color',[0.7,0.7,0.7]);hold on
    
end

for i_xylabel = 1:length(xylabel)
    
    text(i_xylabel*1.2-0.7,-0.5,xylabel{i_xylabel},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    text(-0.8,i_xylabel*1.2-0.6,xylabel{i_xylabel},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',7)
    
end
% text(3,7.5,ligand_all{i_ligand},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'fontsize',9)
xlim([-2,length(xylabel)+2])
ylim([-2,length(xylabel)+2])

legend(h_vec,dataclass,'Location','northeast','box','off');

% medium size figure
paperpos = [0,0,100,100]*5;
papersize = [100,100]*5;
draw_pos = [20,20,60,60]*5;
set(gcf, 'PaperUnits','points')
set(gcf, 'PaperPosition', paperpos,'PaperSize', papersize,'Position',draw_pos)
fig_save_path = '../SubFigures2023/';

set(gca, 'fontsize',7,'XTick',{},'YTick',{})
saveas(gcf,strcat(fig_save_path,fig_save_name),'epsc')
close()

end

