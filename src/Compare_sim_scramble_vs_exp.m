%%
% compare sim vs exp, and scramble sim vs exp
clear all

%% load data
% sim, metrics can be used
load('../raw_data2023/Sim9_all_Comb_ligands_codon_metric.mat')
metrics_sim = metrics;
data_sim = data;

% scrambled sim, metrics can be used
load('../raw_data2023/Sim11_all_comb_scramble_codon_metric.mat')
metrics_scrambled_sim = metrics;
data_scrambled_sim = data;

% exp,metrics can be used
load('Supriya_data_metrics.mat')
metrics_exp = metrics;
data_exp = data;

index_exp_ligand_vects = {[8];%6
    [5];%7
    [9];%8
    [3];%9
    [25];%10
    [19];%11
    [26];%12
    [13,16,18];%13
    [];%14
    [14];%15
    [];%16
    [];%17
    [];%18
    [20];%19
    [];%20
    [21];%21
    [];%22
    [];%23
    [];%24
    [];%25
    [];%26
    [];%27
    [];%28
    [];%29
    [];%30
    []};%31

%% cal codon and sim vs exp diff
codon_list = {'Speed','PeakAmplitude','Duration' ,'TotalActivity', 'EarlyVsLate','OscVsNonOsc'  };
field_list = fieldnames(metrics_exp{1});
% collect_feature_vects_all_minmax_scaled = cell(0);
% metrics_all_minmax_scaled = cell(0);
% collect_feature_vects_all = cell(0);
% metrics_all = cell(0);
% W_diff_sim = cell(0);

for i_cond = 1:length(index_exp_ligand_vects)
    index_exp_ligand = index_exp_ligand_vects{i_cond};
    sti_vec{i_cond} = data_sim.info_ligand(i_cond);
    
    for i_rpc = 1:length(index_exp_ligand)
        clear metrics_cal data_info
        for i_field = 1:length(field_list)
            metrics_cal{1}.(field_list{i_field}) = metrics_exp{index_exp_ligand(i_rpc)}.(field_list{i_field});
            metrics_cal{2}.(field_list{i_field}) = metrics_sim{i_cond}.(field_list{i_field})(1:9:end,:);
            metrics_cal{3}.(field_list{i_field}) = metrics_scrambled_sim{i_cond}.(field_list{i_field})(1:9:end,:);
            
        end
        data_info.info_ligand =[data_exp.info_ligand(index_exp_ligand(i_rpc)),...
            data_sim.info_ligand(i_cond),...
            data_scrambled_sim.info_ligand(i_cond)];
        
        data_info.info_dose_str =[data_exp.info_dose_str(index_exp_ligand(i_rpc)),...
            data_sim.info_dose_str(i_cond),...
            data_scrambled_sim.info_dose_str(i_cond)];
        
        data_info.data_label = {'exp','sim','scramble_sim'};
        
        [collect_feature_vects_all_minmax_scaled{i_cond,i_rpc},metrics_all_minmax_scaled{i_cond,i_rpc}] = calculate_codon_from_metric2023(data_info,metrics_cal); %,  parameter
        [collect_feature_vects_all{i_cond,i_rpc},metrics_all{i_cond,i_rpc}] = calculate_codon_from_metric2023_07_nonminmaxscaled(data_info,metrics_cal); %,  parameter
        
        for i_codon = 1:length(codon_list)
            
            %% exp vs sim
            W_diff_sim{i_cond}(i_codon,i_rpc) = cal_w_dist(collect_feature_vects_all{i_cond,i_rpc}.(codon_list{i_codon}){1},...
                collect_feature_vects_all{i_cond,i_rpc}.(codon_list{i_codon}){2});
            
            
            %% exp vs scramble sim
            W_diff_scrambled_sim{i_cond}(i_codon,i_rpc) = ...
                cal_w_dist(collect_feature_vects_all{i_cond,i_rpc}.(codon_list{i_codon}){1},...
                collect_feature_vects_all{i_cond,i_rpc}.(codon_list{i_codon}){3});
            
            %% minmax rescaled exp vs sim
            W_diff_sim_minmax{i_cond}(i_codon,i_rpc) = cal_w_dist(collect_feature_vects_all_minmax_scaled{i_cond,i_rpc}.(codon_list{i_codon}){1},...
                collect_feature_vects_all_minmax_scaled{i_cond,i_rpc}.(codon_list{i_codon}){2});
            
            %% exp vs scramble sim
            W_diff_scrambled_sim_minmax{i_cond}(i_codon,i_rpc) = ...
                cal_w_dist(collect_feature_vects_all_minmax_scaled{i_cond,i_rpc}.(codon_list{i_codon}){1},...
                collect_feature_vects_all_minmax_scaled{i_cond,i_rpc}.(codon_list{i_codon}){3});
            
        end
        
    end
    % good_fit_pop_wdis(i_ligand,i_dose) = mean(W_diff{i_cond}(:,i_sti));
    
end

W_diff_scrambled_sim_mat = zeros(6,1);
W_diff_sim_mat = zeros(6,1);

for i_cond =1:length(W_diff_scrambled_sim)
    if isempty(W_diff_scrambled_sim{i_cond})
        W_diff_scrambled_sim_mat(:,i_cond) = NaN;
        W_diff_sim_mat(:,i_cond) = NaN;
        W_diff_sim_minmax_mat(:,i_cond) = NaN;
        W_diff_scrambled_sim_minmax_mat(:,i_cond) = NaN;
        
    else
        W_diff_scrambled_sim_mat(:,i_cond) = mean(W_diff_scrambled_sim{i_cond},2);
        W_diff_sim_mat(:,i_cond) = mean(W_diff_sim{i_cond},2);
        W_diff_sim_minmax_mat(:,i_cond) = mean(W_diff_sim_minmax{i_cond},2);
        W_diff_scrambled_sim_minmax_mat(:,i_cond) = mean(W_diff_scrambled_sim_minmax{i_cond},2);
        
    end
end

index = [1:8,10,14,16];
compare_sim = W_diff_scrambled_sim_mat - W_diff_sim_mat;
sum(sign(compare_sim))
sum(sign(compare_sim(:,index)),2)

scramble_sim_good_fit = W_diff_scrambled_sim_mat<=1;
sim_good_fit = W_diff_sim_mat<=1;

index = [1:8,10,14,16];
compare_sim_minmax = W_diff_scrambled_sim_minmax_mat - W_diff_sim_minmax_mat;
sum(sign(compare_sim_minmax))
sum(sign(compare_sim_minmax(:,index)),2)

scramble_sim_good_fit_minmax = W_diff_scrambled_sim_minmax_mat<=0.5;
sim_good_fit_minmax = W_diff_sim_minmax_mat<=0.5;

%% plot
if 1 % W-distance, run this for figure
    clear sti_vec
    fig_save_path = '../SubFigures2023/';
    for i_cond = 1:16
        sti_vec{i_cond} = data_sim.info_ligand(i_cond);
    end
    mymap = [ (0:0.05:1)',(0:0.05:1)',ones(21,1);
        ones(20,1),(0.95:-0.05:0)',(0.95:-0.05:0)'];
    
    figure(2)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,W_diff_sim_minmax_mat,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_wdis_sim_vs_exp'),'epsc')
    close()
    
    figure(1)
    set(gcf, 'PaperUnits','points')
    
    paper_pos = [0,0,450,200];
    paper_size = [450,200];
    set(gcf, 'PaperPosition', paper_pos,'PaperSize',paper_size )%,'Position',draw_pos [20,20,280,280]
    
    h = heatmap( sti_vec,codon_list,W_diff_scrambled_sim_minmax_mat,'Colormap',mymap,'CellLabelColor','none');%'none'
    caxis([0,0.5])
    
    saveas(gcf,strcat(fig_save_path,'codon_diff_wdis_scramble_sim_vs_exp'),'epsc')
    close()
end


%% sub functions
function [w_dis] = cal_w_dist(exp_data,sim_data)

pts = linspace(min(min(exp_data),min(sim_data)),max(max(exp_data),max(sim_data)),50);

[~,~,bw_exp] = ksdensity(exp_data,pts,...
    'Function','pdf');
[~,~,bw_sim] = ksdensity(sim_data,pts,...
    'Function','pdf');%
bw = min(bw_exp,bw_sim);

[f_exp,xi_exp] = ksdensity(exp_data,pts,...
    'Function','pdf','Bandwidth',bw);%,'Bandwidth',bw
[f_sim,xi_sim] = ksdensity(sim_data,pts,...
    'Function','pdf','Bandwidth',bw);%

w_dis = w_distance(exp_data, sim_data, 2);
end
