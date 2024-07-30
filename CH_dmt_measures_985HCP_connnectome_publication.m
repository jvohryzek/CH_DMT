%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2004 Selen Atasoy
% Author: Selen Atasoy (selenatasoy@gmail.com)
% Co-author: Jakub Vohryzek (jakub.vohryzek@upf.edu)
% Date: 21/03/2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
close all
version    = 'v1';
save_now   = 0;

%% folders & files
inputFile       = '/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics';
folder.figures  = fullfile(inputFile,'DMT','Results', 'figures');
folder.dmt      = fullfile('/Users/jakub/Datasets/Psychedelics/DMT/CH_Projections/');
folder.CH       = '/Users/jakub/Datasets/Connectomes';
file.CH         = 'CH_18715_985_HCP_groupconnectome.mat';

file.dmt        = 'DMT_Projections_985HCP_DOT.mat'; % 985fold

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
    Sdist.(mainNames_dmt{f}) = surfFMRI_computeCHspectralDensity(abs(dmt.Projections.(mainNames_dmt{f})), CHBins);
    
    PT.(mainNames_dmt{f})    = dmtFMRI_EnergyOverTimeCH(abs(dmt.Projections.(mainNames_dmt{f})), V);    
    % fitting spectra
end
clear f
% normalise projections by the positive max of the baseline condition
% should it not be absolute value as the negative values might be above -1

normProjectionsBaseline_dmt= dmtFMRI_normalizeProjections(dmt.Projections,'power','DMT_pre_norm'); % gives back normalised projections in relation to your baseline
%% Repertoire Entropy
[h1] = repertoire_entropy_dmt(Sdist)
    

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

%% POWER AND ENERGY

[h4] = dmtFMRI_barCHPower(dmt.Projections,statDist); 
[h5] = dmtFMRI_barCHPower(Eind,statDist);title('Total energy', 'FontSize', 18);

%% save the results
if save_now
   pathFig = '/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_paper/Figures/'

   % Energy
    saveas(h1,strcat(pathFig,'DMT_CH_Repertoire_Entropy_985HCP.png'));
    saveas(h2,strcat(pathFig,'DMT_CH_EnergyBinned_',statDist,'_985HCP.png'));
    saveas(h2extra,strcat(pathFig,'DMT_CH_EnergyBinned_post_',statDist,'_985HCP.png'));
    saveas(h3,strcat(pathFig,'DMT_CH_EnergyDiffBinned_',statDist,'_985HCP.png'));
    saveas(h4,strcat(pathFig,'DMT_CH_PowerBar_',statDist,'_985HCP.png'));
    saveas(h5,strcat(pathFig,'DMT_CH_EnergyBar_',statDist,'_985HCP.png'));

end;   