proj_num_vec = {18,25,26,27,28,[18,26],[18,25],[26,27]};
proj_ligand_vec = {{'TNF'},{'LPS'},{'CpG'},{'polyIC'},{'Pam3CSK'} ,{'TNF','CpG'},{'TNF','LPS'},{'CpG','polyIC'}};
proj_dose_str_vec = {{'10ng/mL'},{'10ng/mL'},{'100nM'},{'100ug/mL'},{'100ng/mL'},{'10ng/mL','100nM'},{'10ng/mL','10ng/mL'},{'100nM','100ug/mL'}};%
proj_dose_val_vec = {{10},{10},{100},{100000},{100},{10,100},{10,10},{100,100000}};

for i_proj = 1:length(proj_ligand_vec)
    
    save_filename = 'para_sampled_NFkBinit';
    
    for i_ligand = 1:length(proj_ligand_vec{i_proj})
        save_filename = strcat(save_filename,'_',proj_ligand_vec{i_proj}{i_ligand},...
            '_',replace(replace(proj_dose_str_vec{i_proj}{i_ligand},'/',''),'.','p'));
    end
    load(strcat(data_save_file_path,save_filename,'.mat'));
    
    a = 1;
end