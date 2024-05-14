function [params, species] = nfkbInitialize_CpGPolyIC_compete()
% This file is automatically generated by updateModel.m from the NFkB whole-cell spreadsheet
% Model URL: multi-stimulus NFkB.xlsx


% PARAMETERS
params(1,1) = 1; % induced IKK activation (general)
params(2,1) = 2e-05; % basal IKK activation
params(3,1) = 2; % IKK inhibition (IKK cycling)
params(4,1) = 18; % IKK renewal (IKK cycling)
params(5,1) = 5e-07; % basal IkBa mRNA synthesis
params(6,1) = 6e-05; % induced IkBa mRNA synthesis
params(6,2) = 2.938; % Hill coefficient for mRNA syn
params(6,3) = 0.1775; % EC50 for mRNA syn
% to delete
params(6,4) = 14; % mRNA transcription/processing/maturation delay
params(7,1) = 0.0577623; % IkBa mRNA degradation
params(8,1) = 30; % IkBa mRNA translation
% to delete
params(8,2) = 1; % Translation/folding delay
params(9,1) = 0.0225; % nuclear import of IkBa
params(9,2) = 3.5; % Volume scale: cytoplasmic/nuclear volume
params(10,1) = 0.6; % nuclear import of NFkB
params(10,2) = params(9,2); % Volume scale: cytoplasmic/nuclear volume
params(11,1) = 0.1575; % nuclear export of IkBa
params(11,2) = 1/params(9,2); % Volume scale: nuclear/cytoplasmic volume
params(12,1) = 0.042; % nuclear export of NFkB
params(12,2) = 1/params(9,2); % Volume scale: nuclear/cytoplasmic volume
params(13,1) = 0; % nuclear import of IkBa-NFkB
%???
params(13,2) = params(9,2); % Volume scale: cytoplasmic/nuclear volume
params(14,1) = 0.828; % nuclear export of IkBa-NFkB
% ???
params(14,2) = 1/params(9,2); % Volume scale: cytoplasmic/nuclear volume
params(15,1) = 0.0770164; % degradation of IkBa
params(16,1) = 0.0770164; % degradation of IkBa (nuc)
params(17,1) = 200; % IkBa-NFkB association
params(18,1) = 200; % IkBa-NFkB association (nuc)
params(19,1) = 0.008; % IkBa-NFkB dissociation
params(20,1) = 0.008; % IkBa-NFkB dissociation (nuc)
params(21,1) = 190; % IKK-IkBa-NFkB association
params(22,1) = 190; % IKK-IkBa association
params(23,1) = 38; % IKK-IkBa-NFkB dissociation
params(24,1) = 38; % IKK-IkBa dissociation
params(25,1) = 2; % Phosphorylation/degradation of complexed IkBa
params(26,1) = 2; % Phosphorylation/degradation of IkBa
params(27,1) = 8.75; % CD14-LPS association
params(27,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(28,1) = 0.07; % CD14-LPS dissociation
params(28,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(29,1) = 0.00112; % CD14 synthesis
params(30,1) = 0.000878; % CD14 degradation
params(31,1) = 5.543; % Assocation of LPS and CD14 in plasma membrane
params(32,1) = 0.0277; % Disassociation of TLR4LPS in the plamsa membrane
params(33,1) = 5.543; % Assocation of LPS and CD14 in the endosome
params(34,1) = 0.0277; % Disassociation of TLR4LPS in the endosome
params(35,1) = 5.25e-05; % Synthesis rate of TLR4
params(36,1) = 0.065681; % Induced endocytosis of CD14
params(37,1) = 0.04; % Recycling of CD14
params(38,1) = 0.028028; % Constitutive endocytosis of TLR4
params(39,1) = 0.5; % Recycling of TLR4
params(40,1) = 0.065681; % Induced endocytosis of activated CD14-TLR4
params(41,1) = 0.04; % Recycling of CD14-TLR4
params(42,1) = 0.07; % Degradation of CD14-LPS
params(43,1) = 0.0653; % Degradation of TLR4
params(44,1) = 0.012; % Degradation of activated TLR4-LPS
params(45,1) = 150; % Activation of MyD88
params(45,2) = 3; % Hill coefficient for MyD88 activation
params(45,3) = 0.012448; % EC50 for MyD88 activation
params(46,1) = 2600; % Deactivation of MyD88
params(47,1) = 6; % Activation of TRIF
params(48,1) = 0.2; % Deactivation of TRIF
params(49,1) = 30; % Activation of TRAF6 by MyD88
params(50,1) = 0.4; % Activation of TRAF6 by TRIF
params(51,1) = 0.125; % Deactivation of TRAF6
params(52,1) = 0.105; % Activation of TAK1 by TRAF6
params(52,2) = 1;
params(53,1) = 0.0116; % TNF degradation
params(54,1) = 8.224e-06; % TNFR synthesis
params(55,1) = 0.02384; % TNFR degradation
params(56,1) = 1100; % Capture of TNF
params(56,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(57,1) = 0.021; % Release of TNF
params(57,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(58,1) = 0.125; % Internalization/degradation of complexed TNFR
params(59,1) = 34.08; % Association of TRAF2/RIP1 with receptor trimer
params(60,1) = 0.03812; % Dissociation of TRAF2/RIP1 with receptor trimer
params(61,1) = 0.125; % Internalization/degradation of complexed TNFR
params(62,1) = 1.875; % Activation of complexed TNFR
params(63,1) = 320.3; % Inactivation of complexed TNFR
params(64,1) = 0.125; % Internalization/degradation of complexed TNFR
params(65,1) = 1889; % Activation of TAK1 by C1
params(65,2) = 1;
params(66,1) = 0.5188; % Inactivation of TAK1
params(67,1) = 18.79; % Activation of IKK by TAK1
params(67,2) = 2; % Hill coefficient for IKK activation
params(67,3) = 0.001116; % EC50 for mRNA syn
params(68,1) = 1e-06; % TLR1/2 synthesis
params(69,1) = 0.0004; % TLR1/2 degradation
params(70,1) = 1; % Association of CD14 and lipoprotein
params(70,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(71,1) = 1.8; % Dissociation of CD14/lipoprotein
params(71,2) = 0.001; % Scale for external (media) vs. internal (cellular) volumes
params(72,1) = 5.543; % NaN
params(73,1) = 0.07; % Degradation of CD14-P3CSK
params(74,1) = 0.02; % NaN
% changed value!!!!
params(75,1) = 0.004; % Degradation of ligand/receptor
params(76,1) = 150; % Activation of MyD88
params(76,2) = 3; % Hill coefficient for MyD88 activation
params(76,3) = 0.0032; % EC50 for MyD88 activation
params(77,1) = 3e-06; % TLR3 synthesis
params(78,1) = 0.0007; % TLR3 degradation
params(79,1) = 0.04; % Poly(I:C) internalization
params(79,2) = 0.03; % EC50 for poly(I:C) internalization
params(79,3) = 1; % Hill coefficient for poly(I:C) internalization
params(79,4) = 0.001; % NaN
params(79,5) = 0.025;%0.025 % 0.8
params(79,6) = 1;
params(79,7) = 0.025;%0.025 % 0.8
params(80,1) = 0.04; % NaN
params(80,2) = 0.001; % Poly(I:C) release
params(81,1) = 0.5; % NaN
params(82,1) = 0.025; % NaN
params(83,1) = 0.0007; % Bound poly(I:C)-TLR3 degradation
params(84,1) = 20; %  Activation of TRIF by bound TLR3
params(85,1) = 2e-06; % TLR9 synthesis
params(86,1) = 0.0004; % TLR9 degradation
params(87,1) = 0.0004; % TLR9 degradation (N terminus fragment)
params(88,1) = 0.015; % CpG internalization
params(88,2) = 0.5; % EC50 for CpG internalization
params(88,3) = 1; % Hill coefficient for CpG internalization
params(88,4) = 0.001; % NaN
params(88,5) = 0.025;%0.025 %3
params(88,6) = 1;
params(88,7) = 0.025;%0.025 %3
params(89,1) = 0.028; % CpG exchange from endosome
params(89,2) = 0.001; % NaN
params(90,1) = 3; % Ligand-receptor association
params(91,1) = 0.5; % Ligand-receptor dissociation
params(92,1) = 3; % Mediated degrdadation of bound TLR9
params(93,1) = 0.0016; % Bound CpG-TLR9 degradation
params(94,1) = 200; % Activation of MyD88
params(94,2) = 3; % Hill coefficient for MyD88 activation
params(94,3) = 0.0032; % EC50 for MyD88 activation
params(95,1) = 0; % degradation of NF-kB (used as proxy for irreversible sequestration)
params(96,1) = 0.33; % IkBat to IkBat_cas1 transformation rate
params(97,1) = 0.33; % IkBat_cas1 to IkBat_cas2 transformation rate
params(98,1) = 0.0577623; % IkBat_cas2 mRNA degradation
params(99,1) = 0.4; %production rate from  NFkBn to NFkBn_cas1
params(100,1) = 0.4; %transformation rate from NFkBn_cas1 to NFkBn_cas2
params(101,1) = 0.4; %degradation rate of NFkBn_cas2
% 0.25??

% SPECIES
species.NAMES = {...
'stim' ...1
'IkBa' ...2
'IkBan' ...3
'IkBaNFkB' ...4
'IkBaNFkBn' ...5
'IkBat' ...6
'IKKIkBaNFkB' ...7
'IKKIkBa' ...8
'NFkB' ...9
'NFkBn' ...10
'IKK_off' ...11
'IKK' ...12
'IKK_i' ...13
'LPS' ...14
'CD14' ...15
'CD14LPS' ...16
'CD14LPSen' ...17
'TLR4' ...18
'TLR4en' ...19
'TLR4LPS' ...20
'TLR4LPSen' ...21
'MyD88_off' ...22
'MyD88' ...23
'TRIF_off' ...24
'TRIF' ...25
'TRAF6_off' ...26
'TRAF6' ...27
'TNF' ...28
'TNFR' ...29
'TNFR_TNF' ...30
'TTR' ...31
'C1_off' ...32
'C1' ...33
'TAK1_off' ...34
'TAK1' ...35
'Pam3CSK' ...36
'TLR2' ...37
'CD14_P3CSK' ...38
'TLR2_P3CSK' ...39
'polyIC' ...40
'polyIC_en' ...41
'TLR3' ...42
'TLR3_polyIC' ...43
'CpG' ...44
'CpG_en' ...45
'TLR9' ...46
'TLR9_CpG' ...47
'TLR9_N' ...48
'IkBat_cas1' ...49
'IkBat_cas2' ...50
'NFkBn_cas1' ...51
'NFkBn_cas2' ...52
};