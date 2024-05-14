% step 1: poly IC parameter 79 traj and osc calculation;
if 0
    parameter_fold_vec = [10.^linspace(-1,1,21)];
    
    Line_wid = 0.75 * ones(length(parameter_fold_vec),1);
    Line_wid(floor(length(parameter_fold_vec)/2)+1) = 1.5;
    color_mapping = [ linspace(1,0.5,floor(length(parameter_fold_vec)/2))',zeros(floor(length(parameter_fold_vec)/2),2);
        0,0,0;
        zeros(floor(length(parameter_fold_vec)/2),2),linspace(0.5,1,floor(length(parameter_fold_vec)/2))'];
    
    data_info.save_file_path = data_save_file_path;
    data_info.species_outputname = {'nucNFkB';'TNFR';'IKK'};
    data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'IKK'}}; % must
    data_info.flag = '';
    data_info.type = 'wt';
    
    sim_info.fold_change_vec = parameter_fold_vec;% 10.^linspace(-1,1,21)];
    sim_info.ligand = {'polyIC'};
    sim_info.dose_str = {'10ug/mL','33ug/mL','100ug/mL'};
    sim_info.dose_val = {10000,33000,100000};
    sim_info.parameter_name = {'params79'};
    para_val = 0.04;
    sim_info.parameter_value_vec = diag(para_val) * sim_info.fold_change_vec;
    
    % sim_data_tbl = NFkB_signaling_para_value_sim(sim_info,data_info);
    
    high_dose = sim_info.dose_val{end};
    
    index_nuc_NFkB = sim_data_tbl.dose_val == high_dose ...
        & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'wt' ...
        & sim_data_tbl.parameter_name == sim_info.parameter_name{1};
    
    i_dose_names =1;
    % NFkBn + IkBaNFkBn
    data.model_sim{i_dose_names} =sim_data_tbl.trajectory(index_nuc_NFkB,1:5:end);
    data.info_ligand{i_dose_names} = sim_data_tbl.ligand(index_nuc_NFkB(1));
    data.info_dose_index{i_dose_names} = 3;
    data.info_dose_str{i_dose_names} = sim_data_tbl.dose_str(index_nuc_NFkB(1));
    data.info_num_cells{i_dose_names} = size(data.model_sim{i_dose_names},1);
    
    data.order{i_dose_names} = (1:data.info_num_cells{i_dose_names})';
    
    data.exp = data.model_sim;
    
    vis_data_field = {'model_sim'};%,'sample'};
    data_label = {'simulation'};%,'sample'};
    
    % [collect_feature_vects,metrics] = calculate_codon(data,vis_data_field,data_label);%,  parameter
    
    % step 2: check whether the order of data/metric/codon is consistent with
    % the order in the sim_data_tbl [checked, the same]
    
    figure(1)
    
    for i_fold_vec = 1:length(sim_info.fold_change_vec)
        index_nuc_NFkB_ind = index_nuc_NFkB   ...
            & sim_data_tbl.parameter_fold_change == sim_info.fold_change_vec(i_fold_vec);
        nucNFkB_traj = sim_data_tbl.trajectory(index_nuc_NFkB_ind,:);
        
        sim_info.fold_change_vec(i_fold_vec)
        metrics{1}.oscpower(i_fold_vec)
        subplot(1,2,1)
        %     plot(1:length(nucNFkB_traj),nucNFkB_traj,'LineWidth',Line_wid(i_fold_vec),'Color',color_mapping(i_fold_vec,:));hold on
        
        subplot(1,2,2)
        nucNFkB_traj2 = data.model_sim{i_dose_names}(i_fold_vec,:);
        %     plot(1:length(nucNFkB_traj2),nucNFkB_traj2,'LineWidth',Line_wid(i_fold_vec),'Color',color_mapping(i_fold_vec,:));hold on
        
    end
    
    % step 3: check osc metrics
    % FramesPerHour = 12;
    sig_stats =get_sig_stats_v202204(data.model_sim{i_dose_names}, 'FS', FramesPerHour);
    
    % step 3.1 check fft: what does this mean:
    % check 5 period; check 6 period; check damped osc; check when freq and
    % amp change the same time
    
    freq_range = [0.33,1];
    
    % freq_period_plot = [0.1,0.2:0.2:1];
    % T_period_plot = 1./ freq_period_plot;
    % figure(2)
    % for i_T_period_plot = 1:length(T_period_plot)
    %         t = 0:5:480;
    %     t = t/60;
    %     time_series = - cos(2*pi*t/T_period_plot(i_T_period_plot)) + 2;%1.5238 1.2308
    % subplot(1,length(T_period_plot),i_T_period_plot)
    %     plot(t, time_series,'Linewidth',1);hold on
    %     title(strcat('freq = ',num2str(freq_period_plot(i_T_period_plot))));
    %
    % xlabel('time')
    % xlim([0,8])
    % ylabel('signal strength')
    % end
    
    
    T_period = 1:0.01:10;
    clear oscpower_T_period
    for i_T_period = 1:length(T_period)
        Fs = 12;  % Sampling frequency
        T_sampling = 1/Fs;  % Sampling period
        Fs = 12;
        t = 0:5:480;
        t = t/60;
        
        time_series = - cos(2*pi*t/T_period(i_T_period)) + 2;%1.5238 1.2308
        
        % fft oscpower calculation
        time_series_lengthen = [time_series, zeros(size(time_series,1),0)];%65536
        L = size(time_series_lengthen,2); % Length of signal
        Y = fft(time_series_lengthen');
        Y = Y';
        P2 = abs(Y/L);
        P1 = P2(:,1:L/2+1);
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        f = Fs*(0:(L/2))/L;
        df = f(2)-f(1);
        sig_stats.oscpower = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
        oscpower_T_period(i_T_period) = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
        
        % Ade's oscpower calculation
        % sig_stats=struct;
        %
        % Data = time_series_lengthen(:, 1:end);
        % Data = fillmissing(Data,'linear', 2, 'EndValues','extrap');
        % smoothData = zeros(size(Data));
        % for i = 1:size(smoothData, 1)
        %     smoothData(i, :) = smooth(Data(i, :), "sgolay");
        % end
        %
        % [psd,fq] = pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
        % sig_stats.fq = fq';
        % sig_stats.psd = transpose(psd./sum(psd,1));
        %
        % psd = transpose(sig_stats.psd) ; fq = transpose(sig_stats.fq);%
        % bp = bandpower(psd,fq,freq_range, 'psd')';
        % sig_stats.oscpower = bp;
        % oscpower_T_period(i_T_period) = sig_stats.oscpower;
        % figure(2)
        % plot(t,time_series);hold on
        % title({'oscpower = ',num2str(sig_stats.oscpower)})
        % xlabel('time (hours)');
    end
    
    figure(3)
    plot(1./T_period,oscpower_T_period/oscpower_T_period(1)); hold on
    xlabel('Signal frequency')
    ylabel('Normalized Oscpower')
    
end
% step 4: check how osc metrics changes for real trajs

if 1
    % step 4-1: pull out the individual traj
    data_save_file_path = '../raw_data/';%_fay_parameter/';
    fig_save_path = '../SubFigures/';%_fay_parameter/';
    
    load(strcat(data_save_file_path,'All_ligand_codon_2.mat'))% data,collect_feature_vects,metrics));
    
    fig_opt.paper_opt.paperpos=[0,0,150,180]*3;
    fig_opt.paper_opt.papersize=[150 180]*3;
    
    % data_proj_num_vec = data_proj_nums{i_ligand};
    %% codon distribution
    % violin_plot_codon_seperate(collect_feature_vects,fig_save_path) % draw each
    % codon and sti seperately, and save
    fig_opt.save_file = strcat(fig_save_path,'All_ligand_Ade_codon_2_distrib');
    field_names_codon = fieldnames(collect_feature_vects);
    index_high_dose = [5,15,25,31,37];
    
    for i_ligand=1:5
        for i_field = 1:length(field_names_codon)
            collect_feature_vects_high_dose.(field_names_codon{i_field}){2*i_ligand-1} = collect_feature_vects.(field_names_codon{i_field}){index_high_dose(i_ligand)};
            collect_feature_vects_high_dose.(field_names_codon{i_field}){2*i_ligand} = collect_feature_vects.(field_names_codon{i_field}){index_high_dose(i_ligand)+1};
        end
    end
    
    cell_index = [50,100,150,200,250];
    
    time_series = metrics{5}.time_series(cell_index,:);%1.5238 1.2308
    time_series = [time_series;metrics{6}.time_series(cell_index,:)];%1.5238 1.2308
    time_series = [time_series;metrics{15}.time_series(cell_index,:)];%1.5238 1.2308
    time_series = [time_series;metrics{16}.time_series(cell_index,:)];%1.5238 1.2308
    
    %     rescale_exp = 1./max(time_series([1:5,11:15],:),[],2)*0.25;
    %     time_series([1:5,11:15],:) = diag(rescale_exp)*time_series([1:5,11:15],:);
    %     time_series([6:10,16:20],:) = diag(rescale_exp)*time_series([6:10,16:20],:);
    
    %     clear time_series
    %     time_series = metrics{5}.time_series(150,:);%1.5238 1.2308
    %     time_series =[time_series; metrics{6}.time_series(150,:)];%1.5238 1.2308
    %
    %     time_series = [time_series;time_series(1:2,:)*2];%1.5238 1.2308
    %     time_series = [time_series;time_series(1:2,:)*3];%1.5238 1.2308
    %     time_series = [time_series;time_series(1:2,:)*4];%1.5238 1.2308
    %     time_series = [time_series;time_series(1:2,:)*5];%1.5238 1.2308
    
    
    %     rescale_exp = 1./max(time_series([1:2:10],:),[],2)*0.25;
    %     time_series_osc([1:2:10],:) = diag(rescale_exp)*time_series([1:2:10],:);
    %     time_series_osc([2:2:10],:) = diag(rescale_exp)*time_series([2:2:10],:);
    
    sti_vec = [3,8,13,16,19];
    time_series = []
    
    for i_sti = sti_vec
        RMSD_vec = sum((data.pred_mode_filter_nan{i_sti} - data.exp_mode_filter_nan{i_sti}).^2,2)/length(data.pred_mode{i_sti});
        nRMSD_vec = RMSD_vec./mean(data.exp_mode_filter_nan{i_sti},2);
        [col1,Ind] = sort(RMSD_vec, 'ascend');
        cdf = cumsum(col1); % Compute cdf
        cdf = cdf/cdf(end); % Normalize
        pct = [0.2,0.5,0.8];
        cell_index = zeros(1,length(pct));
        
        for i_pct =1:length(pct)
            index = find(cdf >= pct(i_pct), 1, 'first');
            cell_index(i_pct) = Ind(index);
        end
        
        time_series = [time_series; data.exp_mode_filter_nan{i_sti}(cell_index,:)];%1.5238 1.2308
        time_series = [time_series;data.pred_mode_filter_nan{i_sti}(cell_index,:)];%1.5238 1.2308
    end
    
    %     exp_index =[1:length(cell_index),2*length(cell_index)+1:3*length(cell_index)];
    %     model_index = [length(cell_index)+1:2*length(cell_index),3*length(cell_index)+1:4*length(cell_index)];
    %     rescale_exp = 1./max(time_series(exp_index,:),[],2)*0.25;
    %     time_series_osc(exp_index,:) = diag(rescale_exp)*time_series(exp_index,:);
    %     time_series_osc(model_index,:) = diag(rescale_exp)*time_series(model_index,:);
    
    
    
    
    % fft oscpower calculation
    time_series_lengthen = [time_series, zeros(size(time_series,1),0)];%65536
    %         time_series_lengthen = [time_series_osc, zeros(size(time_series_osc,1),0)];%65536
    
    %% oscpower fft
    if 1
        rescale_exp = 1./max(time_series(:,:),[],2)*0.25;
        time_series_osc = diag(rescale_exp)*time_series(:,:);
        
        time_series_lengthen = [time_series_osc, zeros(size(time_series_osc,1),0)];%65536
        Fs = 12;  % Sampling frequency
        T_sampling = 1/Fs;  % Sampling period
        freq_range = [0.33,1];
        Fs = 12;
        L = size(time_series_lengthen,2); % Length of signal
        Y = fft(time_series_lengthen');
        Y = Y';
        P2 = abs(Y/L);
        P1 = P2(:,1:L/2+1);
        P1(:,2:end-1) = 2*P1(:,2:end-1);
        f = Fs*(0:(L/2))/L;
        df = f(2)-f(1);
        sig_stats.oscpower = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
        oscpower_real_traj = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
    end
    
    
    %% oscpwoer welch
    % Ade's oscpower calculation
    if 0
        sig_stats=struct;
        Fs = 12;  % Sampling frequency
        freq_range = [0.33,1];
        
        Data = time_series_lengthen(:, 1:end);
        Data = fillmissing(Data,'linear', 2, 'EndValues','extrap');
        smoothData = zeros(size(Data));
        for i = 1:size(smoothData, 1)
            smoothData(i, :) = smooth(Data(i, :), "sgolay");
        end
        
        [psd,fq] = pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
        sig_stats.fq = fq';
        sig_stats.psd = transpose(psd./sum(psd,1));
        
        psd = transpose(sig_stats.psd) ; fq = transpose(sig_stats.fq);%
        bp = bandpower(psd,fq,freq_range, 'psd')';
        sig_stats.oscpower = bp;
        oscpower_real_traj = bp;
    end
    
    %     for i_cell = 1:size(time_series,1)/2
    %         subplot(1,size(time_series,1)/2, i_cell)
    %         cell_exp_index = i_cell*2-1;
    %         cell_model_index = i_cell*2;
    %         plot((0:size(time_series,2)-1)/12,time_series(cell_exp_index,:),'k','LineWidth',1);hold on
    %         plot((0:size(time_series,2)-1)/12,time_series(cell_model_index,:),'r','LineWidth',1);hold on
    %         legend({strcat('exp-',num2str(oscpower_real_traj(cell_exp_index))),...
    %             strcat('model-',num2str(oscpower_real_traj(cell_model_index)))})
    %
    %         ylim([0,0.25])
    %         xlim([0,8])
    %     end
    
    title_str = {'TNF-10ng/mL','LPS-100ng/mL'};
    title_str = {'TNF-10ng/mL','LPS-100ng/mL','CpG-1ug/mL','pIC-100ug/mL','P3C4-100ng/mL',};
    figure(1)
    close
    figure(1)
    ligand_vec = [1,2,3,4,5];
    for i_ligand = 1:length(ligand_vec)
        for i_cell = 1:length(cell_index)
            
            cell_exp_index = i_cell+(i_ligand-1)*2*length(cell_index);
            cell_model_index = i_cell+(i_ligand-1)*2*length(cell_index) + length(cell_index);
            subplot(length(ligand_vec),length(cell_index),(i_ligand-1)*length(cell_index) + i_cell)
            plot((0:size(time_series,2)-1)/12,time_series(cell_exp_index,:),'k','LineWidth',1);hold on
            plot((0:size(time_series,2)-1)/12,time_series(cell_model_index,:),'r','LineWidth',1);hold on
            legend({strcat('exp-',num2str(oscpower_real_traj(cell_exp_index))),...
                strcat('model-',num2str(oscpower_real_traj(cell_model_index)))})
            legend('boxoff')
            
            ylim([0,0.25])
            xlim([0,8])
            title(title_str{i_ligand})
        end
    end
    
    
    if 0
        for i_sti = 1:length(metrics)
            time_series_lengthen = metrics{i_sti}.time_series;
            
            if 0
                sig_stats=struct;
                Fs = 12;  % Sampling frequency
                freq_range = [0.33,1];
                
                Data = time_series_lengthen(:, 1:end);
                Data = fillmissing(Data,'linear', 2, 'EndValues','extrap');
                smoothData = zeros(size(Data));
                for i = 1:size(smoothData, 1)
                    smoothData(i, :) = smooth(Data(i, :), "sgolay");
                end
                
                [psd,fq] = pwelch(smoothData',size(Data,2),10,256,Fs,'one-sided','psd');
                sig_stats.fq = fq';
                sig_stats.psd = transpose(psd./sum(psd,1));
                
                psd = transpose(sig_stats.psd) ; fq = transpose(sig_stats.fq);%
                bp = bandpower(psd,fq,freq_range, 'psd')';
                sig_stats.oscpower = bp;
                metrics{i_sti}.oscpower = bp;
                vects0{i_sti} = bp;
            end
            
            if 1
                Fs = 12;  % Sampling frequency
                T_sampling = 1/Fs;  % Sampling period
                freq_range = [0.33,1];
                Fs = 12;
                rescale_exp = 1./max(time_series_lengthen,[],2)*0.25;
                time_series_lengthen = diag(rescale_exp)*time_series_lengthen;
                
                L = size(time_series_lengthen,2); % Length of signal
                Y = fft(time_series_lengthen');
                Y = Y';
                P2 = abs(Y/L);
                P1 = P2(:,1:L/2+1);
                P1(:,2:end-1) = 2*P1(:,2:end-1);
                f = Fs*(0:(L/2))/L;
                df = f(2)-f(1);
                metrics{i_sti}.oscpower = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
                vects0{i_sti} = sum(P1(:,f>=freq_range(1) & f<=freq_range(2)).^2*df,2);
            end
        end
        
        all_data = cell2mat(vects0(:)); %%% normalize each codeword
        all_data = (all_data-prctile(all_data,0.5))/(prctile(all_data,99.5)-prctile(all_data,0.5)); %normalized except those extrema
        count = 1;
        for k = 1:length(vects)
            vects{k} = all_data(count:(size(vects{k}, 1)+count-1), :);
            count = count + size(vects{k}, 1);
        end
        collect_feature_vects.OscVsNonOsc = vects';
        
        fig_opt.save_file = strcat(fig_save_path,'All_ligand_normalizedosc_codon_2_distrib');
        violin_plot_codon(collect_feature_vects,fig_opt) % draw all codon and sti in one
        
        
        field_names_codon = fieldnames(collect_feature_vects);
        index_high_dose = [5,15,25,31,37];
        
        for i_ligand=1:5
            for i_field = 1:length(field_names_codon)
                collect_feature_vects_high_dose.(field_names_codon{i_field}){2*i_ligand-1} = collect_feature_vects.(field_names_codon{i_field}){index_high_dose(i_ligand)};
                collect_feature_vects_high_dose.(field_names_codon{i_field}){2*i_ligand} = collect_feature_vects.(field_names_codon{i_field}){index_high_dose(i_ligand)+1};
            end
        end
        
        fig_opt.save_file = strcat(fig_save_path,'high_ligand_normalizedosc_codon_2_distrib');
        violin_plot_codon(collect_feature_vects_high_dose,fig_opt) % draw all codon and sti in one
    end
end