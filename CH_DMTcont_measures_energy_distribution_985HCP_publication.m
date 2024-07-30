%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jakub Vohryzek (jakub.vohryzek@upf.edu)
% Original script - author: Selen Atasoy (selenatasoy@gmail.com)
% Calculating the time-resolved energy spectrum
% Date: 14/03/2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_now   = 0;

%% folders & files
inputFile       = '/Users/jakub/Datasets/Psychedelics/DMT_cont/CH_Projections';
inputFile_CH    = '/Users/jakub/Datasets/Connectomes';
file_CH         = 'CH_18715_985_HCP_groupconnectome'; % original CH
file_dmt        = 'DMTcont_Projections_985HCP_DOT.mat'; % original projections
%% load connectome harmonics
load(fullfile(inputFile_CH, file_CH));
CH = Connectome_Harmonics;
V  = Eigenvalues(1:10000,1); % taking the first 100000 values
%% load dmt data
dmt = load(fullfile(inputFile, file_dmt));
mainNames_dmt = fieldnames(dmt.Projections);

%% define the dimensions
numSbj = size(dmt.Projections.(mainNames_dmt{1}), 1); % dmt number of subjects
numTp  = size(dmt.Projections.(mainNames_dmt{1}), 2); % dmt number of timepoints
numCH  = size(dmt.Projections.(mainNames_dmt{1}), 3); % number of harmonics
for f = 1:length(mainNames_dmt);
    dmt.Projections.(mainNames_dmt{f})     = dmt.Projections.(mainNames_dmt{f})(:,:,1:10000);
end
%% compute the energy and power distributions
% calculating Edist for 15 bins across CH space
numBins = 15;
[CHBins,CHBinCenters] = surfFMRI_binCHinds(CH(1,1:numCH-8715), numBins, 1); % 1 for logspace, 0 for linspace
numBins = length(CHBinCenters);

for f=1:length(mainNames_dmt)
    Edist.(mainNames_dmt{f}) = surfFMRI_computeCHenergyDistribution(abs(dmt.Projections.(mainNames_dmt{f})), V, CHBins); % dim 15x237x15 (sbj,tps,bins)
end
clear f
%% labels and colors
mainNames = fieldnames(Edist);
mainNamesLabels = ({'DMT PRE';'DMT POST';'PCB PRE';'PCB POST'});
mainNamesLabelsDiff = ({'DMT POST - DMT PRE';'PCB POST - PCB PRE'});
cMap = dmtFMRI_colorMap(size(mainNames,1));

%% intensity ratings
load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_cont/intensityRatings_DMTcont.mat')
figure,plot(mean(intensityRatings_dmt_cont,1)','LineWidth',2)
hold on
plot(mean(intensityRatings_pcb_cont)','LineWidth',2)
ylabel('Intensity Rating');xlabel('Minutes');xlim([1 29])
legend({'DMT','PLACEBO'})

%% avearging for pre and post and without normalising to reflect the changes in DMT non-continuous
dt = 30;

% overall global pattern
% figure,bar(squeeze(mean(mean(Edist.DMT_all_norm)))-squeeze(mean(mean(Edist.PCB_all_norm))))

% pre-post pattern
figure,
subplot(2,2,1);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,240:480,:))))-squeeze(mean(mean(Edist.DMT_all_norm(:,1:240,:)))));ylim([-0.03 0.03])
title('DMTpost vs. DMTpre');
subplot(2,2,2);bar(squeeze(mean(mean(Edist.PCB_all_norm(:,240:480,:))))-squeeze(mean(mean(Edist.PCB_all_norm(:,1:240,:)))));ylim([-0.03 0.03])
title('PCBpost vs. PCBpre');
%pre-cond and post-cond
subplot(2,2,3);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,1:240,:))))-squeeze(mean(mean(Edist.PCB_all_norm(:,1:240,:)))));ylim([-0.03 0.03])
title('DMTpre vs. PCBpre');
subplot(2,2,4);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,240:480,:))))-squeeze(mean(mean(Edist.PCB_all_norm(:,240:480,:)))));ylim([-0.03 0.03])
title('DMTpost vs. PCBpost');

%% the mean at pre and post for both groups
DMT240mean = squeeze(mean(mean(mean(Edist.DMT_all_norm(:,1:240,:)))));
DMT480mean = squeeze(mean(mean(mean(Edist.DMT_all_norm(:,240:480,:)))));
PCB240mean = squeeze(mean(mean(mean(Edist.PCB_all_norm(:,1:240,:)))));
PCB480mean = squeeze(mean(mean(mean(Edist.PCB_all_norm(:,240:480,:)))));

figure,bar([DMT240mean DMT480mean PCB240mean PCB480mean])
%% avearging for pre and post and normalising to reflect the changes in DMT non-continuous
figure,
subplot(2,2,1);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,240:480,:))))./DMT480mean-squeeze(mean(mean(Edist.DMT_all_norm(:,1:240,:))))./DMT240mean);ylim([-0.3 0.3])
title('DMTpost vs. DMTpre');
subplot(2,2,2);bar(squeeze(mean(mean(Edist.PCB_all_norm(:,240:480,:))))./PCB480mean-squeeze(mean(mean(Edist.PCB_all_norm(:,1:240,:))))./PCB240mean);ylim([-0.3 0.3])
title('PCBpost vs. PCBpre');
%pre-cond and post-cond
subplot(2,2,3);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,1:240,:))))./DMT240mean-squeeze(mean(mean(Edist.PCB_all_norm(:,1:240,:))))./PCB240mean);ylim([-0.3 0.3])
title('DMTpre vs. PCBpre');
subplot(2,2,4);bar(squeeze(mean(mean(Edist.DMT_all_norm(:,240:480,:))))./DMT480mean-squeeze(mean(mean(Edist.PCB_all_norm(:,240:480,:))))./PCB480mean);ylim([-0.3 0.3])
title('DMTpost vs. PCBpost');
%% normalising measures by the group mean at every timepoint
for i =1:28
    % the mean at every timepoint
    DMTcont_mean(i) = squeeze(mean(mean(mean(Edist.DMT_all_norm(:,dt*(i-1)+1:dt*i,:)))))
    PCBcont_mean(i) = squeeze(mean(mean(mean(Edist.PCB_all_norm(:,dt*(i-1)+1:dt*i,:)))))
    % normalised measures by the mean at every timepoint
    Edist_DMT_min(i,:,:) = squeeze(mean(Edist.DMT_all_norm(:,dt*(i-1)+1:dt*i,:),2))./DMTcont_mean(i);
    Edist_PCB_min(i,:,:) = squeeze(mean(Edist.PCB_all_norm(:,dt*(i-1)+1:dt*i,:),2))./PCBcont_mean(i);
    Edist_diff_normalised(i,:,:) = squeeze(Edist_DMT_min(i,:,:))-squeeze(Edist_PCB_min(i,:,:));
end
figure
subplot(121);bar(DMTcont_mean)
subplot(122);bar(PCBcont_mean)
%% energy difference spectrum per minute average across subjects
figure
numMin = 28

for t=1:numMin % first 28 minutes
    subplot(4,7,t)
    bar(squeeze(mean(Edist_diff_normalised(t,:,:),2)));axis square;axis off
    ylim([-0.2 0.2])
    title(['Minute ',num2str(t)])

end
% %% by subject
% figure
% % 
% for sbj=1:14
%     for t=1:numMin % first 28 minutes
%         subplot(4,7,t)
%         bar(squeeze(mean(Edist_diff_normalised(t,sbj,:),2)))
%         ylim([-0.2 0.2])
% 
%     end
%     pause
% end
%% %%%%%%%%%%%%%%%%%%%% spectrum difference %%%%%%%%%%%%%%%%%%%%
% ground truth
gt = load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT/chPowerEnergy_DMT.mat')

% ground truth spectrum
Edist_diff_ground_truth = squeeze(mean(mean(gt.Edist.DMT_post1_norm-gt.Edist.DMT_pre_norm,2),1));
%%
for i=1:14
    for t=1:numMin
        meas_normalised(i,t) = corr(squeeze(Edist_diff_normalised(t,i,:)),Edist_diff_ground_truth,'Type', 'Spearman');
    end
end

%%
for s = 1:14
    [h_spect p_spect] = corr(intensityRatings_dmt_cont(s,:)'-intensityRatings_pcb_cont(s,:)',meas_normalised(s,:)','type','Spearman');
    h_meas_DMT(s) = h_spect;
    p_meas_DMT(s) = p_spect;

end

%%
figure
yyaxis left
plot(1:1:28,mean(intensityRatings_dmt_cont-intensityRatings_pcb_cont,1)','LineWidth',2)
hold on
ylabel('Intensity Rating')
yyaxis right
%plot(squeeze(nanmean(meas_min,2))','m-','LineWidth',2)
plot(squeeze(nanmean(meas_normalised,1))','m-','LineWidth',2)

plot([8 8],[-0.3 0.4],'r-','LineWidth',2)
xlim([0 28]);%ylim([-0.15 0.2])
set(gca,"XTick",0:1:28,'XTickLabel',0:1:28,'ycolor','m')
legend('DMT - PCB intensity rating','DMT - PCB Power spectrum','injection','Location','southeast')
ylabel('Power Spectrum Difference');xlabel('Minutes')


%% plotting
Sbj_names={'Sbj1','Sbj2','Sbj3','Sbj4','Sbj5','Sbj6','Sbj7',...
    'Sbj8','Sbj9','Sbj10','Sbj11','Sbj12','Sbj13','Sbj14'};
figure
subplot(4,1,[1 2 3])
%for net=1:7
    barh(1:14,h_meas_DMT,'k')%,'FaceColor',color_vecs(net,:))
    hold on
%end
set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
title('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
xlim([-0.6 0.8])

box off
for net=1:14
    if p_meas_DMT(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(h_meas_DMT(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(h_meas_DMT(1,net)),' Pval: ',num2str(p_meas_DMT(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end
subplot(4,1,4)

boxplot(h_meas_DMT,'Group Average')
camroll(270)
ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(h_meas_DMT,0)

%% %%%% LZ and entropy and spectrum diff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% visualising the LZ regressors
% load('RegressorLZInterpscrubbed.mat')
% figure,
% plot(regDMT','k');hold on;
% plot(regPCB','r');hold on;
% plot(mean(regPCB,1)','y','LineWidth',2);hold on;
% plot(mean(regDMT,1)','b','LineWidth',2)
% 
% load('RegressorLZInterpscrubbedConvolvedAvg.mat')
% figure,
% plot(RegDMT2','k');hold on;
% plot(RegPCB2','r');hold on;
% plot(mean(RegPCB2,1)','y','LineWidth',2);hold on;
% plot(mean(RegDMT2,1)','b','LineWidth',2)
% 
% figure,
% plot(regdiff','k');hold on;
% plot(mean(regdiff,1)','y','LineWidth',2);hold on;
% plot(Regdiff','k');hold on;
% plot(mean(Regdiff,1)','r','LineWidth',2);hold on;
% 
% %% reduce to 28 minutes bins
% %%
% dt = 30
% for i=1:28
%     regDMT_min(i,:) = squeeze(mean(regDMT(:,dt*(i-1)+1:dt*i),2))';
%     RegDMT2_min(i,:) = squeeze(mean(RegDMT2(:,dt*(i-1)+1:dt*i),2))';
%     regPCB_min(i,:) = squeeze(mean(regPCB(:,dt*(i-1)+1:dt*i),2))';
%     RegPCB2_min(i,:) = squeeze(mean(RegPCB2(:,dt*(i-1)+1:dt*i),2))';
%     regdiff_min(i,:) = squeeze(mean(regdiff(:,dt*(i-1)+1:dt*i),2))';
%     Regdiff_min(i,:) = squeeze(mean(Regdiff(:,dt*(i-1)+1:dt*i),2))';
% end
% %% visualising the LZ regressors
% 
% figure,
% plot(regdiff_min,'k');hold on;
% plot(mean(regdiff_min,2),'y','LineWidth',2);hold on;
% plot(Regdiff_min,'k');hold on;
% plot(mean(Regdiff_min,2),'r','LineWidth',2);hold on;
% %%
% for s = 1:14
%     [h_spect p_spect] = corrcoef(regDMT_min(:,s)',Edist_diff_normalised(:,s)');
%     h_LZvsSpectdiff_DMT(s) = h_spect(1,2);
%     p_LZvsSpectdiff_DMT(s) = p_spect(1,2);
% 
% end
% %% plotting
% 
% figure
% subplot(4,1,[1 2 3])
% 
% barh(1:14,h_LZvsSpectdiff_DMT,'k')
% hold on
% 
% set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
% set(gca,'YTickLabelRotation',45)
% xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
% %ylim([-.4 1.1])
% box off
% for net=1:14
%     if p_LZvsSpectdiff_DMT(1,net)<(0.05)
%         % To plot an asterisk when significant
%         plot(h_LZvsSpectdiff_DMT(1,net)*1.1,net,'*k')
%             string2annot = {strcat('Corr: ',num2str(h_LZvsSpectdiff_DMT(1,net)),' Pval: ',num2str(p_LZvsSpectdiff_DMT(1,net)))};
%         %text(5,0.4,string2annot,'FontSize',20)
%     end
% end
% 
% set(gca,'YGrid','on','FontSize',14)
% subplot(4,1,4)
% 
% boxplot(h_LZvsSpectdiff_DMT,'Group Average')
% camroll(270)
% %ylim([-0.3 0.7])
% set(gca,'YTickLabel',[]);
% 
% [h p] = ttest(h_LZvsSpectdiff_DMT,0)
% 
% 
% %% correlations
% 
% figure
% for s=1:14
% subplot(3,5,s)
% x3 = regDMT_min(:,s); % different data but the same mean5HT2A_aalsymm(:,1); % taking the symmetric version
% y3 = Edist_diff_normalised(:,s);
% 
% % correlation
% [RHO,PVAL] = corr(x3,y3,'Type','Spearman');
% 
% p = polyfit(x3,y3,1);
% % Evaluate the fitted polynomial p and plot:
% f = polyval(p,x3);
% plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
% xlabel({'Intensity Ratings'})
% ylabel('Entropy')
% title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
% set(gca, 'FontSize',14)
% 
% end
