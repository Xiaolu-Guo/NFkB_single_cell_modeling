function [] =  draw_traj_heatmap_co_sti_CpG_polyIC_competetation(data_file,fig_save_path,data_save_file_path_1)
% For Guo et al. Figure 4F, heatmap CpG-polyIC stimulation with competition 
% 
% tested 05/12/2024, Matlab 2020a

%% NFkB ordered by NFkB total activity
if 1
    clear data
    
    load(strcat(data_save_file_path_1,data_file))% data,collect_feature_vects,metrics));
    
    %[~,data.order{1}] = sort(max(data.model_sim{1}(1:9:end,:),[],2),'descend');
    [~,data.order{1}] = sort(sum(data.model_sim{1}(1:9:end,:),2),'descend');

    data_fields = {'model_sim','exp'};
    for i_data_fields = 1:length(data_fields)
        data.(data_fields{i_data_fields}){1} = data.(data_fields{i_data_fields}){1}(1:9:end,:);
    end
    
    vis_data_field = {'model_sim'};
    fig_opts.fig_name = 'Debug_NFkB_ordered_NFkB_';
    fig_opts.fig_savepath = fig_save_path;
    fig_opts.color_limits = [-0.01,0.25];
    
    max(max(data.model_sim{1}))
    plot_traj_heatmap(data,vis_data_field,fig_opts)
end


