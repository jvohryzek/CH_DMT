%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2004 Selen Atasoy
% Author: Selen Atasoy (selenatasoy@gmail.com)
% Co-author: Jakub Vohryzek (jakub.vohryzek@queens.ox.ac.uk)
% Date: 03/07/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_now   = 1; % '0' or '1'
psiloData  = 0; % '0' or '1'
lsdData    = 0; % '0' or '1'
excludeSbj = 'All'; % 'Most_Stringent', 'All','DMTcont','Least_Stringent' subjects with missing brain regions
pathFig = ['/Users/jakub/Figures/Collaboration_Kringelbach/Connectome_Harmonics/27_10_2021_excSbj_',excludeSbj,'/'];

%% folders & files
inputFile       = '/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics';
%folder.figures  = fullfile(inputFile,'DMT','Results', 'figures');
folder.dmt      = fullfile('/Users/jakub/Datasets/Psychedelics/DMT/CH_Projections');
folder.psilo    = fullfile(inputFile,'PSILO','psiloProjectionCH');
folder.lsd      = fullfile(inputFile,'LSD','Rest1'); % the first Rest1 of the lsd condition 
folder.CH       = fullfile('/Users/jakub/Datasets/Connectomes');

file.CH         = 'CH_sym_average_weighted.mat';
file.dmt        = 'DMT_Projections_sym_prob_DOT.mat';
file.psilo      = 'PSILO_Projections_sym_average_weighted_DOT.mat';
file.lsd        = 'LSD_Projections_sym_average_weighted_DOT.mat';

%% load the files  
% load connectome harmonics
load(fullfile(folder.CH, file.CH));
% dmt data
switch excludeSbj
    case 'All'
        dmt = load(fullfile(folder.dmt, file.dmt));
   case 'DMTcont'
        sbj2keep = [2 3 4 5 7 8 9 10 11 12 13 14 16 17]
        % sbj2exclude = [2 3 9 10]
        dmt = load(fullfile(folder.dmt, file.dmt));
        dmt.Projections.DMT_pre_norm = dmt.Projections.DMT_pre_norm(sbj2keep,:,:);
        dmt.Projections.DMT_post1_norm = dmt.Projections.DMT_post1_norm(sbj2keep,:,:);
        dmt.Projections.PCB_pre_norm = dmt.Projections.PCB_pre_norm(sbj2keep,:,:);
        dmt.Projections.PCB_post1_norm = dmt.Projections.PCB_post1_norm(sbj2keep,:,:);

    case 'Most_Stringent'
        sbj2keep = [1 8 11 13 14 15 16 17]
        % sbj2exclude = [2 3 4 5 6 7 9 10 12]
        dmt = load(fullfile(folder.dmt, file.dmt));
        dmt.Projections.DMT_pre_norm = dmt.Projections.DMT_pre_norm(sbj2keep,:,:);
        dmt.Projections.DMT_post1_norm = dmt.Projections.DMT_post1_norm(sbj2keep,:,:);
        dmt.Projections.PCB_pre_norm = dmt.Projections.PCB_pre_norm(sbj2keep,:,:);
        dmt.Projections.PCB_post1_norm = dmt.Projections.PCB_post1_norm(sbj2keep,:,:);
    case 'Least_Stringent'
        sbj2keep = [1 4 5 6 7 8 11 12 13 14 15 16 17]
        % sbj2exclude = [2 3 9 10]
        dmt = load(fullfile(folder.dmt, file.dmt));
        dmt.Projections.DMT_pre_norm = dmt.Projections.DMT_pre_norm(sbj2keep,:,:);
        dmt.Projections.DMT_post1_norm = dmt.Projections.DMT_post1_norm(sbj2keep,:,:);
        dmt.Projections.PCB_pre_norm = dmt.Projections.PCB_pre_norm(sbj2keep,:,:);
        dmt.Projections.PCB_post1_norm = dmt.Projections.PCB_post1_norm(sbj2keep,:,:);

end
mainNames_dmt = fieldnames(dmt.Projections);
% define the dimensions
numSbj = size(dmt.Projections.(mainNames_dmt{1}), 1); % dmt number of subjects
numTp  = size(dmt.Projections.(mainNames_dmt{1}), 2); % dmt number of timepoints
numCH  = size(dmt.Projections.(mainNames_dmt{1}), 3); % number of harmonics

% psilo data
if psiloData
    psilo = load(fullfile(folder.psilo, file.psilo));
    % version part of psilodep
    dmt.Projections.placebo_psilo = psilo.Projections.PCB_norm;
    dmt.Projections.psilo = psilo.Projections.PSILO_norm;

    % version on its own
    % psilo.Projections.psilo = psilo.Projections.PSILO_norm;
    % psilo.Projections.placebo = psilo.Projections.PCB_norm;
    % psilo.Projections = rmfield(psilo.Projections,'PCB_norm');
    % psilo.Projections = rmfield(psilo.Projections,'PSILO_norm');
end

% lsd data
if lsdData
    lsd = load(fullfile(folder.lsd, file.lsd));
    % version part of psilodep
    dmt.Projections.placebo_lsd = lsd.Projections.PCB_norm;
    dmt.Projections.lsd = lsd.Projections.LSD_norm;

    % version on its own
    % psilo.Projections.psilo = psilo.Projections.PSILO_norm;
    % psilo.Projections.placebo = psilo.Projections.PCB_norm;
    % psilo.Projections = rmfield(psilo.Projections,'PCB_norm');
    % psilo.Projections = rmfield(psilo.Projections,'PSILO_norm');
end

% defining the main fields
mainNames = fieldnames(dmt.Projections);
if psiloData
    numSbj_psilo = size(dmt.Projections.(mainNames{5}), 1); % psilo number of subjects
    numTp_psilo = size(dmt.Projections.(mainNames{5}), 2); % psilo number of timepoints
end
if lsdData
    numSbj_lsd = size(dmt.Projections.(mainNames{7}), 1); % lsd number of subjects
    numTp_lsd = size(dmt.Projections.(mainNames{7}), 2); % lsd number of timepoints
end
%% compute the energy and power distributions
% calculating log scale for 15 bins across CH space
numBins = 15;
[CHBins,CHBinCenters] = surfFMRI_binCHinds(CH(1,1:numCH), numBins, 1); % 1 for logspace, 0 for linspace
numBins = length(CHBinCenters);

for f=1:length(mainNames)
    Edist.(mainNames{f}) = surfFMRI_computeCHenergyDistribution(abs(dmt.Projections.(mainNames{f})), V, CHBins); % dim 15x237x15 (sbj,tps,bins)
    Eind.(mainNames{f})  = surfFMRI_computeCHenergyIndividual(abs(dmt.Projections.(mainNames{f})), V);
    Sdist.(mainNames{f}) = surfFMRI_computeCHspectralDensity(abs(dmt.Projections.(mainNames{f})), CHBins);
    
    PT.(mainNames{f})    = dmtFMRI_EnergyOverTimeCH(abs(dmt.Projections.(mainNames{f})), V);    
    % fitting spectra
    Efit.(mainNames{f})  = dmtFMRI_fitSemiLogEnergy(Edist.(mainNames{f}), CHBinCenters, 'mean');
end
clear f
% normalise projections by the positive max of the baseline condition
% should it not be absolute value as the negative values might be above -1

normProjectionsBaseline_dmt = dmtFMRI_normalizeProjections(dmt.Projections,'power','DMT_pre_norm'); % gives back normalised projections in relation to your baseline
if psiloData
    normProjectionsBaseline_psilo = dmtFMRI_normalizeProjections(dmt.Projections,'power','pcb_psilo'); % gives back normalised projections in relation to your baseline
end
if lsdData
    normProjectionsBaseline_lsd = dmtFMRI_normalizeProjections(dmt.Projections,'power','pcb_lsd'); % gives back normalised projections in relation to your baseline
end
%% Repertoire Entropy
% TODO: so far in this script implementation of repetori entropy for all
% conditions only (DMT,PCB,LSD,PSILO)
if psiloData && lsdData
    mainNamesLabelsDiff = ({'DIFF_DMT';'DIFF_PCB';'DIFF_PSILO';'DIFF_LSD'});

    [ent_norm,entsum_norm] = repertoire_Entropy(Sdist,mainNames);
    for f = 1:length(mainNames);
        ent_norm_mean.(mainNames{f})     = mean(ent_norm.(mainNames{f}),2);
        entsum_norm_mean.(mainNames{f})  = mean(entsum_norm.(mainNames{f}),2);
    end
    idx=0;
    for f = 1:2:length(mainNames)
        idx=idx+1;
        ent_norm_mean_diff.(mainNamesLabelsDiff{idx})     = mean(ent_norm.(mainNames{f+1}),2) - mean(ent_norm.(mainNames{f}),2);
        entsum_norm_mean_diff.(mainNamesLabelsDiff{idx})  = mean(entsum_norm.(mainNames{f+1}),2) - mean(entsum_norm.(mainNames{f}),2);
    end
    [h10] = repertoire_entropy_plot(entsum_norm_mean,entsum_norm_mean_diff,mainNames)
else
    [h10,entStats] = repertoire_entropy_dmt(Sdist)
end
%% correlation with subjective scores
sbjScores = readtable('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_paper/subjectives_scores_publication.xlsx')
richExp = sbjScores.How_rich_WasYourConciousExperience;
idSbj = [4,6,7,8,9,10,11,12,13,14,15,16,17]; % asubjects with non-nan values in the scores
[x p] = corr(squeeze(entStats(2,idSbj)'-entStats(1,idSbj)'),richExp(idSbj),'Type','Spearman')
figure,plot(entStats(2,idSbj)-entStats(1,idSbj),richExp(idSbj),'*','LineWidth',4);ylabel('Richness of experience');xlabel('CH Repertoire Entropy')
title(['Correlation = ', num2str(x(1,2)),' p-val = ',num2str(p(1,2))])
%
%% correlation with subjective scores
% sbjScores = readtable('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_paper/subjectives_scores_publication.xlsx')
% richExp = sbjScores.How_rich_WasYourConciousExperience;
% richExp_rmout = rmoutliers(richExp);
% [x p] = corrcoef(squeeze(entStats(2,[6,7,8,10,11,12,13,14,15,16,17])),richExp([6,7,8,10,11,12,13,14,15,16,17])')
% figure,plot(entStats(2,6:end),richExp(6:end),'*','LineWidth',4);ylabel('Richness of experience');xlabel('CH Repertoire Entropy')
% title(['Correlation = ', num2str(x(1,2)),' p-val = ',num2str(p(1,2))])
% figure,plot(entStats(2,6:end)-entStats(1,6:end),richExp(6:end),'*','LineWidth',4);ylabel('Richness of experience');xlabel('CH Repertoire Entropy')

%% FIGURES
%% DISTRIBUTION AND DIFFERENCE IN ENERGY
statDist = 'Subjects'; % statDist = 'SubjectsxTimepoints';

[h2] = dmtFMRI_barCHEnergyDistribution(Edist, CHBinCenters,statDist);
[h3] = dmtFMRI_barCHEnergyDistributionDiff(Edist, CHBinCenters,statDist);

%% plot - checking how temporally consistent the patterns are
%%
numMin = 24;
dt = 10;
for i=1:numMin
    % meas_min(i,:) = squeeze(mean(meas(:,dt*(i-1)+1:dt*i),2))';
    %Edist_diff_min(i,:,:) = squeeze(mean(Edist_diff(:,dt*(i-1)+1:dt*i,:),2));
    Edist_DMTpre_min(i,:,:) = squeeze(mean(Edist.DMT_pre_norm(:,dt*(i-1)+1:dt*i,:),2));
    Edist_PCBpre_min(i,:,:) = squeeze(mean(Edist.PCB_pre_norm(:,dt*(i-1)+1:dt*i,:),2));
    Edist_DMTpost_min(i,:,:) = squeeze(mean(Edist.DMT_post1_norm(:,dt*(i-1)+1:dt*i,:),2));
    Edist_PCBpost_min(i,:,:) = squeeze(mean(Edist.PCB_post1_norm(:,dt*(i-1)+1:dt*i,:),2));
    Edist_DMTdiff_min(i,:,:) = Edist_DMTpost_min(i,:,:) - Edist_DMTpre_min(i,:,:) ;
    Edist_PCBdiff_min(i,:,:) = Edist_PCBpost_min(i,:,:) - Edist_PCBpre_min(i,:,:);
end
%%
figure
for i=1:2
    for t=1:numMin % first 16 minutes
        subplot(2,numMin,t+(i-1)*numMin)
        if i ==1
         bar(squeeze(mean(Edist_DMTdiff_min(t,:,:),2)));axis off;axis square;ylim([-0.02 0.02])
        elseif i==2
            bar(squeeze(mean(Edist_PCBdiff_min(t,:,:),2)));axis off;axis square;ylim([0 0.3])

        elseif i== 3
           bar(squeeze(mean(Edist_PCB_min(t,:,:),2)));axis off;axis square;ylim([0 0.3])
        elseif i== 4
           bar(squeeze(mean(Edist_DMT_min(t,:,:),2))-squeeze(mean(Edist_PCB_min(t,:,:),2)));axis off;axis square;ylim([-0.02 0.02])

        end
    end
end
%% Extra - comparing post and pre conditions
Edist_extra = struct;
Edist_extra.(mainNames{4}) = Edist.(mainNames{4});
Edist_extra.(mainNames{2}) = Edist.(mainNames{2});
Edist_extra.(mainNames{3}) = Edist.(mainNames{3});
Edist_extra.(mainNames{1}) = Edist.(mainNames{1});

[h2extra] = dmtFMRI_barCHEnergyDistribution(Edist_extra, CHBinCenters,statDist);
[h3extra] = dmtFMRI_barCHEnergyDistributionDiff(Edist_extra, CHBinCenters,statDist);
%% 
figure,
subplot(2,2,1);bar(squeeze(mean(mean(Edist.DMT_post1_norm)))-squeeze(mean(mean(Edist.DMT_pre_norm))));ylim([-0.03 0.03])
title('DMTpost vs. DMTpre');
subplot(2,2,2);bar(squeeze(mean(mean(Edist.PCB_post1_norm)))-squeeze(mean(mean(Edist.PCB_pre_norm))));ylim([-0.03 0.03])
title('PCBpost vs. PCBpre');
%pre-cond and post-cond
subplot(2,2,3);bar(squeeze(mean(mean(Edist.DMT_pre_norm)))-squeeze(mean(mean(Edist.PCB_pre_norm))));ylim([-0.03 0.03])
title('DMTpre vs. PCBpre');
subplot(2,2,4);bar(squeeze(mean(mean(Edist.DMT_post1_norm)))-squeeze(mean(mean(Edist.PCB_post1_norm))));ylim([-0.03 0.03])
title('DMTpost vs. PCBpost');
%% POWER AND ENERGY

[h8] = dmtFMRI_barCHPower(dmt.Projections,statDist); 
[h5] = dmtFMRI_barCHPower(Eind,statDist);title('Total energy', 'FontSize', 18);

%% NORMALISED POWER DISTRIBUTION AND ENERGY DISTRIBUTION
if psiloData
    [h9psilo] = dmtFMRI_plotCHProjectionDistribution(normProjectionsBaseline_psilo, 1:18715);
end
if lsdData
    [h9lsd] = dmtFMRI_plotCHProjectionDistribution(normProjectionsBaseline_lsd, 1:18715);
end
[h9dmt] = dmtFMRI_plotCHProjectionDistribution(normProjectionsBaseline_dmt, 1:18715);
[h7] = dmtFMRI_plotCHEnergyDistribution(PT);

%% log-normal dist

[h12] = dmtFMRI_plotCHEnergySemiLog(Edist, CHBinCenters, 'mean'); % same as [h3] mean values
[h13] = dmtFMRI_plotCHEnergySemiLog(Edist, CHBinCenters,'std'); % same as [h3] stddev values

%% save the results
if save_now
   %pathFig = [folder.figures,'/'];
   if ~exist(pathFig, 'dir')
       mkdir(pathFig)
   end
   % Energy
   if psiloData && lsdData
    saveas(h2,strcat(pathFig,'CH_EnergyBinned_psilolsd_',statDist,'.png'));
    saveas(h3,strcat(pathFig,'CH_EnergyDiffBinned_psilolsd_',statDist,'.png'));

    saveas(h5,strcat(pathFig,'CH_EnergyBar_psilolsd_',statDist,'.png'));
    saveas(h8,strcat(pathFig,'CH_PowerBar_psilolsd_',statDist,'.png'));
    
    saveas(h10,strcat(pathFig,'CH_Repeertoire_entropy_psilolsd_',statDist,'.png'));
   else
    saveas(h2,strcat(pathFig,'CH_EnergyBinned_',statDist,'.png'));
    saveas(h2extra,strcat(pathFig,'CH_EnergyBinned_post_condition',statDist,'.png'));

    saveas(h3,strcat(pathFig,'CH_EnergyDiffBinned_',statDist,'.png'));

    saveas(h5,strcat(pathFig,'CH_EnergyBar_',statDist,'.png'));
    saveas(h8,strcat(pathFig,'CH_PowerBar_',statDist,'.png'));
    
    saveas(h10,strcat(pathFig,'CH_Repeertoire_entropy',statDist,'.png'));

   end
    %%
    if psiloData
        saveas(h9psilo,strcat(pathFig,'CH_ProjectionDistribution_psilo.png'));
    end
    if lsdData
        saveas(h9lsd,strcat(pathFig,'CH_ProjectionDistribution_lsd.png'));
    end
    saveas(h9dmt,strcat(pathFig,'CH_ProjectionDistribution_dmt.png'));

    saveas(h7,strcat(pathFig,'CH_EnergyDistribution.png'));

    saveas(h12,strcat(pathFig,'CH_MeanEnergy.png'));
    saveas(h13,strcat(pathFig,'CH_EnergyFluctuations.png'));
            
            

end;   