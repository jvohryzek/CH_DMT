%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2004 Selen Atasoy
% Author: Selen Atasoy (selenatasoy@gmail.com)
% Co-author: Jakub Vohryzek (jakub.vohryzek@queens.ox.ac.uk)
% Date: 03/07/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
close all
version    = 'v1';
save_now   = 0;

%% folders & files
inputFile       = '/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics';
folder.figures  = fullfile(inputFile,'DMT','Results', 'figures');
folder.dmt      = fullfile('/Users/jakub/Datasets/Psychedelics/DMT/CH_Projections/');
folder.CH       = '/Users/jakub/Datasets/Connectomes'; % fullfile(inputFile);
file.CH         = 'CH_18715_HCP_MGH_32fold_groupconnectome.mat'; % 'CH_sym_average_weighted.mat';
%file.CH         = 'CH_18715_985_HCP_groupconnectome.mat'; % 'CH_sym_average_weighted.mat';

%file.dmt        = 'DMT_Projections_sym_prob_DOT.mat'; % original
file.dmt        = 'DMT_Projections_32fold_DOT.mat'; % 32fold

%% load the files  
% load connectome harmonics
load(fullfile(folder.CH, file.CH));
CH = Connectome_Harmonics;
V  = Eigenvalues(1:10000,1); % taking the first 100000 values
% dmt data

dmt = load(fullfile(folder.dmt, file.dmt));
mainNames_dmt = fieldnames(dmt.Projections);
% define the dimensions
numSbj = size(dmt.Projections.(mainNames_dmt{1}), 1); % dmt number of subjects
numTp  = size(dmt.Projections.(mainNames_dmt{1}), 2); % dmt number of timepoints
numCH  = size(dmt.Projections.(mainNames_dmt{1}), 3); % number of harmonics
for f = 1:length(mainNames_dmt);
    dmt.Projections.(mainNames_dmt{f})     = dmt.Projections.(mainNames_dmt{f})(:,:,1:10000);
end
%% compute the energy and power distributions
% calculating log scale for 15 bins across CH space
numBins = 15;
[CHBins,CHBinCenters] = surfFMRI_binCHinds(CH(1,1:numCH-8715), numBins, 1); % 1 for logspace, 0 for linspace
numBins = length(CHBinCenters);

for f=1:length(mainNames_dmt)
    Edist.(mainNames_dmt{f}) = surfFMRI_computeCHenergyDistribution(abs(dmt.Projections.(mainNames_dmt{f})), V, CHBins); % dim 15x237x15 (sbj,tps,bins)
    Eind.(mainNames_dmt{f})  = surfFMRI_computeCHenergyIndividual(abs(dmt.Projections.(mainNames_dmt{f})), V);
    % Pind.(mainNames{f})  = abs(dmt.Projections.(mainNames{f}));
    Sdist.(mainNames_dmt{f}) = surfFMRI_computeCHspectralDensity(abs(dmt.Projections.(mainNames_dmt{f})), CHBins);
    
    PT.(mainNames_dmt{f})    = dmtFMRI_EnergyOverTimeCH(abs(dmt.Projections.(mainNames_dmt{f})), V);    
    % fitting spectra
    %Efit.(mainNames{f})  = dmtFMRI_fitSemiLogEnergy(Edist.(mainNames{f}), CHBinCenters, 'mean');
end
clear f
% normalise projections by the positive max of the baseline condition
% should it not be absolute value as the negative values might be above -1

normProjectionsBaseline_dmt= dmtFMRI_normalizeProjections(dmt.Projections,'power','DMT_pre_norm'); % gives back normalised projections in relation to your baseline
%% Repertoire Entropy
[h10] = repertoire_entropy_dmt(Sdist)
    

%% FIGURES
%% DISTRIBUTION AND DIFFERENCE IN ENERGY
statDist = 'Subjects'; % statDist = 'SubjectsxTimepoints';

[h2] = dmtFMRI_barCHEnergyDistribution(Edist, CHBinCenters,statDist);

[h3] = dmtFMRI_barCHEnergyDistributionDiff(Edist, CHBinCenters,statDist);
%% Extra - comparing post and pre conditions
Edist_extra = struct;
Edist_extra.(mainNames_dmt{4}) = Edist.(mainNames_dmt{4});
Edist_extra.(mainNames_dmt{2}) = Edist.(mainNames_dmt{2});
Edist_extra.(mainNames_dmt{3}) = Edist.(mainNames_dmt{3});
Edist_extra.(mainNames_dmt{1}) = Edist.(mainNames_dmt{1});

[h2extra] = dmtFMRI_barCHEnergyDistribution(Edist_extra, CHBinCenters,statDist);

% comparing pre_dmt as baseline to the rest
% EdistDepHealthy.psilo = Edist.psilo
% EdistDepHealthy.DMT_pre_norm1 = Edist.DMT_pre_norm
% 
% EdistDepHealthy.DMT_pre_norm2 = Edist.DMT_pre_norm
% EdistDepHealthy.PCB_post1_norm = Edist.PCB_post1_norm
% 
% EdistDepHealthy.DMT_pre_norm3 = Edist.DMT_pre_norm
% EdistDepHealthy.DMT_post1_norm = Edist.DMT_post1_norm
% 
% [h2] = psilodepFMRI_barCHEnergyDistributionDiff(EdistDepHealthy, CHBinCenters,statDist);
% [h3] = psilodepFMRI_barCHEnergyDistribution(EdistDepHealthy, CHBinCenters,statDist);

% OLDv: [h3a] = lsdFMRI_barCHEnergyDistribution(Edist, CHBinCenters);
% OLDv: [h2a] = lsdFMRI_barCHEnergyDistributionDiff(Edist, CHBinCenters);
%% POWER AND ENERGY

[h8] = dmtFMRI_barCHPower(dmt.Projections,statDist); 
[h5] = dmtFMRI_barCHPower(Eind,statDist);title('Total energy', 'FontSize', 18);
% [hadd] = lsdFMRI_showCHProjections(dmt.Projections, [2:51]);
% [hadd] = lsdFMRI_showCHProjections(Eind, [2:51]);

% unnecessary - but a nice check of the 
% data[h5] = dmtFMRI_barCHPower(dmtFMRI_normalizeProjections(dmt.Projections, V, 'energy'));

%% NORMALISED POWER DISTRIBUTION AND ENERGY DISTRIBUTION

[h9dmt] = dmtFMRI_plotCHProjectionDistribution(normProjectionsBaseline_dmt, 1:10000);

[h7] = dmtFMRI_plotCHEnergyDistribution(PT);

%% log-normal dist
% OLDv: [h12] = lsdFMRI_plotCHEnergySemiLog(Edist, CHBinCenters, 'mean'); % same as [h3] mean values
% OLDv: [h13] = lsdFMRI_plotCHEnergySemiLog(Edist, CHBinCenters,'std'); % same as [h3] stddev values

[h12] = dmtFMRI_plotCHEnergySemiLog(Edist, CHBinCenters, 'mean'); % same as [h3] mean values
[h13] = dmtFMRI_plotCHEnergySemiLog(Edist, CHBinCenters,'std'); % same as [h3] stddev values

%% save the results
if save_now
   %pathFig = [folder.figures,'/'];
   pathFig = '/Users/jakub/Figures/Collaboration_Kringelbach/Connectome_Harmonics/21_03_2021/'

   % Energy

    saveas(h2,strcat(pathFig,'CH_EnergyBinned_',statDist,'_32fold.png'));
    %saveas(h2extra,strcat(pathFig,'CH_EnergyBinned_post_condition',statDist,'_32fold.png'));

    saveas(h3,strcat(pathFig,'CH_EnergyDiffBinned_',statDist,'_32fold.png'));

    saveas(h5,strcat(pathFig,'CH_EnergyBar_',statDist,'_32fold.png'));
    saveas(h8,strcat(pathFig,'CH_PowerBar_',statDist,'_32fold.png'));
    
    saveas(h10,strcat(pathFig,'CH_Repertoire_Entropy_32fold.png'));
    %%
    saveas(h9dmt,strcat(pathFig,'CH_ProjectionDistribution_dmt_32fold.png'));

    saveas(h7,strcat(pathFig,'CH_EnergyDistribution_32fold.png'));

    saveas(h12,strcat(pathFig,'CH_MeanEnergy_32fold.png'));
    saveas(h13,strcat(pathFig,'CH_EnergyFluctuations_32fold.png'));
            
            
    %%
    % saveas(h14,strcat(path,'CH_EnergyGaussian_mean.jpg'));       
    % saveas(h15,strcat(path,'CH_EnergyGaussian_std.jpg'));
            
    
end;   