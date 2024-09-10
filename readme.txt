This is instruction for running codes of Guo et. al. "Modeling heterogeneous signaling dynamics of macrophages reveals principles of information transmission in stimulus responses"

Figure S1:
(1) run_me_publish_nc.m
run section "parameter sens analysis" for both run the simulation and visualizing the results. (set values in the if 语句 为1）.

Figure 2:
(1) run_me_publish_nc.m
run section "rescaling the data to SI for monolix software input" to rescale the experimental data to S.I. units and change the data format to monolix software input format. 
(2) run monolix codes in the following order:
(2-a) run XGESN2023001.mlxtran to estimate the parameter for the ODE parameter distribution in the core module. 
(2-b) run XGESN2023002.mlxtran to estimate the parameter for the ODE parameter distribution in the TNF module, with the core module population distribution fixed. 
(2-c) run XGESN2023003.mlxtran to estimate the parameter for the ODE parameter distribution in the LPS module, with the core module population distribution fixed. 
(2-d) run XGESN2023004.mlxtran to estimate the parameter for the ODE parameter distribution in the CpG module, with the core module population distribution fixed. 
(2-e) run XGESN2023005.mlxtran to estimate the parameter for the ODE parameter distribution in the PolyIC module, with the core module population distribution fixed. 
(2-f) run XGESN2023006.mlxtran to estimate the parameter for the ODE parameter distribution in the Pam3CSK module, with the core module population distribution fixed. 
(3) run_me_publish_nc.m
run section "calculate signaling codon" to calculate the signaling codons and save the results.
(4) run_me_publish_nc.m
run section "Heatmaps of traj, W-dist of signaling codon distri." to get Figure 2ABC
(5) to calculate the machine learning classification based on the signaling codons 
(5-a) run_me_publish_nc.m
run section "codon prep for MI calculation & machine learning format" to get the machine learning format signaling codons dataset, which is the input for the classification algorithm.
(5-b) run ??.py to to train the XGBoost (?) Random Forest model and make predictions to get the the confusion matrix.
(5-c) run_me_publish_nc.m
run section "stimulation classificatio" to get the results for machine learning classification results, to get Figure 2D

Figure S2: 
(1) run_me_publish_nc.m
run section "metrics of good fitting, RMSD, signaling codon distri, and W-dist" to get figure S2A-E.
(2) to calculate the machine learning classification based on the signaling codons 
(2-a) run_me_publish_nc.m
run section "codon prep for MI calculation & machine learning format" to get the machine learning format signaling codons dataset, which is the input for the classification algorithm.
(2-b) run ??.py to to train the XGBoost (?) Random Forest model and make predictions to get the the confusion matrix.
(2-c) run_me_publish_nc.m
run section "stimulation classificatio" to get the results for machine learning classification results, to get Figure 2F-G

Figure 3:


