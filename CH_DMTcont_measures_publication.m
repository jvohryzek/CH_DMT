%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2004 Selen Atasoy
% Author: Jakub Vohryzek (jakub.vohryzek@upf.edu)
% Original script - author: Selen Atasoy (selenatasoy@gmail.com)

% Date: 20/11/2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all
version    = 'v1';
save_now   = 0;

%% folders & files
inputFile       = '/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_cont/Results';
inputFile_CH    = '/Users/jakub/Datasets/Connectomes';
file_CH         = 'CH_18715_HCP_MGH_32fold_groupconnectome.mat';
file_dmt        = 'DMT_Projections_sym_prob_DOT.mat'; % original

%% load the files  
% load connectome harmonics
load(fullfile(inputFile_CH, file_CH));
CH =  Connectome_Harmonics;
V  = Eigenvalues(1:10000,1); % taking the first 100000 values
%% load dmt data
dmt = load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_cont/Results/DMT_Projections_sym_prob_DOT.mat', 'Projections')
dmt = load(fullfile(inputFile, file_dmt));
mainNames_dmt = fieldnames(dmt.Projections);

%% define the dimensions
numSbj = size(dmt.Projections.(mainNames_dmt{1}), 1); % dmt number of subjects
numTp  = size(dmt.Projections.(mainNames_dmt{1}), 2); % dmt number of timepoints
numCH  = size(dmt.Projections.(mainNames_dmt{1}), 3); % number of harmonics
for f = 1:length(mainNames_dmt)
    dmt.Projections.(mainNames_dmt{f})     = dmt.Projections.(mainNames_dmt{f})(:,:,1:10000);
end
%% compute the energy and power distributions
% calculating Sdist for 15 bins across CH space
numBins = 15;
[CHBins,CHBinCenters] = surfFMRI_binCHinds(CH(1,1:numCH-8715), numBins, 1); % 1 for logspace, 0 for linspace
numBins = length(CHBinCenters);

%
for f=1:length(mainNames_dmt)
    Edist.(mainNames_dmt{f}) = surfFMRI_computeCHenergyDistribution(abs(dmt.Projections.(mainNames_dmt{f})), V, CHBins); % dim 15x237x15 (sbj,tps,bins)
    Sdist.(mainNames_dmt{f}) = surfFMRI_computeCHspectralDensity(abs(dmt.Projections.(mainNames_dmt{f})), CHBins);
end
clear f
%% intensity ratings
load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_cont/intensityRatings_DMTcont.mat')
figure,plot(mean(intensityRatings_dmt_cont,1)','LineWidth',2)
hold on
plot(mean(intensityRatings_pcb_cont)','LineWidth',2)
ylabel('Intensity Rating');xlabel('Minutes');xlim([1 29])
legend({'DMT','PLACEBO'})
%% Repertoire Entropy
% [h10] = repertoire_entropy_dmt_cont(Sdist)
%%
% labels and colors
mainNames = fieldnames(Sdist);
mainNamesLabels = ({'DMT PRE';'DMT POST';'PCB PRE';'PCB POST'});
mainNamesLabelsDiff = ({'DMT POST - DMT PRE';'PCB POST - PCB PRE'});


cMap =  dmtFMRI_colorMap(size(mainNames,1));


%% entropy calculation
sbjSdist=struct;
for f=1:length(mainNames)
    for sbj=1:size(Sdist.DMT_all_norm,1)
        sbjSdist = squeeze(Sdist.(mainNames{f})(sbj,:,:));
        Prob_PowCH = sbjSdist';
        Prob_PowCHsum = (sbjSdist'./sum(sbjSdist'));

        ent = -sum(Prob_PowCH.*log2(Prob_PowCH));
        entsum = -sum(Prob_PowCHsum.*log2(Prob_PowCHsum));

        ent_norm(f,sbj,:) = ent./log2(size(Prob_PowCH,1));
        entsum_norm(f,sbj,:) = entsum./log2(size(Prob_PowCHsum,1));
        
        clear sbjSdist Prob_PowCH Prob_PowCHsum ent entsum
    end
end
% figure,plot(squeeze(mean(ent_norm,2))');legend(mainNamesLabels)
% figure,violinplot(squeeze(mean(ent_norm,3))',mainNamesLabels);ylabel('Repertoire Entropy')
entStats = squeeze(mean(entsum_norm,3));
%%
% figure
% plot(squeeze(entsum_norm(1,:,:))','r')
% hold on
% plot(squeeze(entsum_norm(2,:,:))','k')
% plot([240 240],[0.84 0.92],'g','LineWidth',2)
%%
figure
yyaxis left
plot(30:30:840,mean(intensityRatings_dmt_cont,1)','LineWidth',2)
hold on
plot(30:30:840,mean(intensityRatings_pcb_cont,1)','LineWidth',2)
ylabel('Intensity Rating')
yyaxis right
plot(squeeze(nanmean(entsum_norm(1,:,:),2))','m-','LineWidth',2)
plot(squeeze(nanmean(entsum_norm(2,:,:),2))','g-','LineWidth',2)
plot([240 240],[0.84 0.91],'r-','LineWidth',2)
xlim([0 840])
set(gca,"XTick",0:30:840,'XTickLabel',0:28,'ycolor','m')
legend('DMT intensity rating','PCB intensity rating','DMT','PCB','injection')
ylabel('Entropy');xlabel('Minutes')

%%
% figure
% plot(squeeze(nanmean(entsum_norm(1,:,1:end),2))','r')
% hold on
% plot(squeeze(nanmean(entsum_norm(2,:,1:end),2))','k')
% hold on
% plot(squeeze(nanmean(entsum_norm(1,:,1:end),2))'-squeeze(nanmean(entsum_norm(2,:,1:end),2))','b')
% plot([240 240],[0.84 0.92],'g','LineWidth',2)
% plot(mean(intensityRatings_pcb_cont)','g','LineWidth',2)
% 
% %hold on
%plot(1:30:840,mean(intensityRatings_dmt_cont,1)'./max(mean(intensityRatings_dmt_cont,1)'))

%%
% tmp_DMT = corrcoef(squeeze(entsum_norm(1,:,1:end))')-eye(14);
% tmp_PCB = corrcoef(squeeze(entsum_norm(2,:,1:end))')-eye(14);
% tmp_DMT_triu = triu(tmp_DMT,1);
% tmp_PCB_triu = triu(tmp_PCB,1);
% 
% figure
% subplot(121)
% imagesc(tmp_DMT,[0 0.2]);axis square;colorbar
% subplot(122)
% imagesc(tmp_PCB,[0 0.2]);axis square;colorbar
% 
% figure
% histogram(tmp_DMT(:),20)
% hold on
% histogram(tmp_PCB(:),20)


%%
dt = 30
for i=1:28
    entsum_norm_dmt_min(i,:) = squeeze(mean(entsum_norm(1,:,dt*(i-1)+1:dt*i),3))';
    entsum_norm_pcb_min(i,:) = squeeze(mean(entsum_norm(2,:,dt*(i-1)+1:dt*i),3))';
end
%%
figure
yyaxis left
plot(1:1:28,mean(intensityRatings_dmt_cont,1)','LineWidth',2)
hold on
plot(1:1:28,mean(intensityRatings_pcb_cont,1)','LineWidth',2)
ylabel('Intensity Rating')
yyaxis right
plot(squeeze(nanmean(entsum_norm_dmt_min,2))','m-','LineWidth',2)
plot(squeeze(nanmean(entsum_norm_pcb_min,2))','g-','LineWidth',2)
plot([8 8],[0.84 0.92],'r-','LineWidth',2)
xlim([0 28]);ylim([0.86 0.92])
set(gca,"XTick",0:1:28,'XTickLabel',0:1:28,'ycolor','m')
legend('DMT intensity rating','PCB intensity rating','DMT','PCB','injection')
ylabel('Entropy');xlabel('Minutes')
%%
for s = 1:14
    [h, p]= corrcoef(intensityRatings_dmt_cont(s,:),entsum_norm_dmt_min(:,s)')
    %tmp_PCB(s) = corrcoef(intensityRatings_pcb_cont(s,:)+0.01,entsum_norm_pcb_min(:,s)')
    corr_DMT_intensity(s) = h(1,2)
    pval_DMT_intensity(s) = p(1,2)

end
%% plotting
Sbj_names={'Sbj1','Sbj2','Sbj3','Sbj4','Sbj5','Sbj6','Sbj7',...
    'Sbj8','Sbj9','Sbj10','Sbj11','Sbj12','Sbj13','Sbj14'};

figure
subplot(4,1,[1 2 3])

barh(1:14,corr_DMT_intensity,'k')%,'FaceColor',color_vecs(net,:))
hold on

set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
title('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
xlim([0 0.9])

box off
for net=1:14
    if pval_DMT_intensity(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(corr_DMT_intensity(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(corr_DMT_intensity(1,net)),' Pval: ',num2str(pval_DMT_intensity(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(corr_DMT_intensity,'Group Average')
camroll(270)
ylim([0 0.9])
set(gca,'YTickLabel',[]);

[h p] = ttest(corr_DMT_intensity,0)
%% correlations

figure
for s=1:14
subplot(3,5,s)
x3 = intensityRatings_dmt_cont(s,:)'; % different data but the samee mean5HT2A_aalsymm(:,1); % taking the symmetric version
y3 = entsum_norm_dmt_min(:,s);

% correlation
[RHO,PVAL] = corr(x3,y3,'Type','Spearman')

p = polyfit(x3,y3,1);
% Evaluate the fitted polynomial p and plot:
f = polyval(p,x3);
plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'Intensity Ratings'})
ylabel('Entropy')
title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
set(gca, 'FontSize',14)

end
%%
% figure
% subplot(121)
% imagesc(tmp_DMT,[0 0.2]);axis square;colorbar
% subplot(122)
% imagesc(tmp_PCB,[0 0.2]);axis square;colorbar


%% %%%%%%%%%%%%%%%%%%%% spectrum difference %%%%%%%%%%%%%%%%%%%%
% ground truth
gt = load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT/chPowerEnergy_DMT.mat')
%%
% Sdist_diff = Sdist.(mainNames_dmt{1})-Sdist.(mainNames_dmt{2});
 Edist_diff = Edist.DMT_all_norm-Edist.PCB_all_norm;

 % ground truth spectrum
Edist_diff_ground_truth = squeeze(mean(mean(gt.Edist.DMT_post1_norm-gt.Edist.DMT_pre_norm,2),1));
Edist_diff_ground_truth;
for i=1:14
    for t=1:840
    meas(i,t) = corr2(squeeze(Edist_diff(i,t,:)),Edist_diff_ground_truth);
    end
end

%% intensity ratings
%load('/Users/jakub/Matlab/Collaboration_Atasoy/ConnectomeHarmonics/DMT_cont/intensityRatings_DMTcont.mat')
%figure,plot(mean(intensityRatings_dmt_cont-intensityRatings_pcb_cont,1)')

%%
dt = 30
for i=1:28
    % meas_min(i,:) = squeeze(mean(meas(:,dt*(i-1)+1:dt*i),2))';
    Edist_diff_min(i,:,:) = squeeze(mean(Edist_diff(:,dt*(i-1)+1:dt*i,:),2));
    Edist_DMT_min(i,:,:) = squeeze(mean(Edist.DMT_all_norm(:,dt*(i-1)+1:dt*i,:),2));
    Edist_PCB_min(i,:,:) = squeeze(mean(Edist.PCB_all_norm(:,dt*(i-1)+1:dt*i,:),2));

end
%%
for i=1:14
    for t=1:28
    meas_min(t,i) = corr2(squeeze(Edist_diff_min(t,i,:)),Edist_diff_ground_truth);
    end
end
%%
for s = 1:14
    [h_spect p_spect] = corrcoef(intensityRatings_dmt_cont(s,:)-intensityRatings_pcb_cont(s,:),meas_min(:,s)');
    h_meas_DMT(s) = h_spect(1,2);
    p_meas_DMT(s) = p_spect(1,2);

end

%%
figure
yyaxis left
plot(1:1:28,mean(intensityRatings_dmt_cont-intensityRatings_pcb_cont,1)','LineWidth',2)
hold on
ylabel('Intensity Rating')
yyaxis right
plot(squeeze(nanmean(meas_min,2))','m-','LineWidth',2)
plot([8 8],[-0.15 0.2],'r-','LineWidth',2)
xlim([0 28]);%ylim([-0.15 0.2])
set(gca,"XTick",0:1:28,'XTickLabel',0:1:28,'ycolor','m')
legend('DMT - PCB intensity rating','DMT - PCB Power spectrum','injection')
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
xlim([-0.3 0.7])

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
%%
figure
yyaxis left
plot(30:30:840,mean(intensityRatings_dmt_cont-intensityRatings_pcb_cont,1)','LineWidth',2)
hold on

ylabel('Intensity Rating')
yyaxis right
plot(squeeze(nanmean(meas,1))','m-','LineWidth',2)
plot([240 240],[-0.2 0.3],'r-','LineWidth',2)
ylim([-0.2 0.3]);xlim([0 840])
set(gca,"XTick",0:30:840,'XTickLabel',0:28,'ycolor','m')
legend('DMT intensity rating','PCB intensity rating','DMT','PCB','injection')
ylabel('Entropy');xlabel('Minutes')

%%%%%% LZ and entropy and spectrum diff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% visualising the LZ regressors
load('RegressorLZInterpscrubbed.mat')
figure,
plot(regDMT','k');hold on;
plot(regPCB','r');hold on;
plot(mean(regPCB,1)','y','LineWidth',2);hold on;
plot(mean(regDMT,1)','b','LineWidth',2)

load('RegressorLZInterpscrubbedConvolvedAvg.mat')
figure,
plot(RegDMT2','k');hold on;
plot(RegPCB2','r');hold on;
plot(mean(RegPCB2,1)','y','LineWidth',2);hold on;
plot(mean(RegDMT2,1)','b','LineWidth',2)

figure,
plot(regdiff','k');hold on;
plot(mean(regdiff,1)','y','LineWidth',2);hold on;
plot(Regdiff','k');hold on;
plot(mean(Regdiff,1)','r','LineWidth',2);hold on;

%% reduce to 28 minutes bins
%%
dt = 30
for i=1:28
    regDMT_min(i,:) = squeeze(mean(regDMT(:,dt*(i-1)+1:dt*i),2))';
    RegDMT2_min(i,:) = squeeze(mean(RegDMT2(:,dt*(i-1)+1:dt*i),2))';
    regPCB_min(i,:) = squeeze(mean(regPCB(:,dt*(i-1)+1:dt*i),2))';
    RegPCB2_min(i,:) = squeeze(mean(RegPCB2(:,dt*(i-1)+1:dt*i),2))';
    regdiff_min(i,:) = squeeze(mean(regdiff(:,dt*(i-1)+1:dt*i),2))';
    Regdiff_min(i,:) = squeeze(mean(Regdiff(:,dt*(i-1)+1:dt*i),2))';
end
%% visualising the LZ regressors

figure,
plot(regdiff_min,'k');hold on;
plot(mean(regdiff_min,2),'y','LineWidth',2);hold on;
plot(Regdiff_min,'k');hold on;
plot(mean(Regdiff_min,2),'r','LineWidth',2);hold on;
%%
for s = 1:14
    [h_spect p_spect] = corrcoef(regDMT_min(:,s)',meas_min(:,s)');
    h_LZvsSpectdiff_DMT(s) = h_spect(1,2);
    p_LZvsSpectdiff_DMT(s) = p_spect(1,2);

end
%% plotting

figure
subplot(4,1,[1 2 3])

barh(1:14,h_LZvsSpectdiff_DMT,'k')
hold on

set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if p_LZvsSpectdiff_DMT(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(h_LZvsSpectdiff_DMT(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(h_LZvsSpectdiff_DMT(1,net)),' Pval: ',num2str(p_LZvsSpectdiff_DMT(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(h_LZvsSpectdiff_DMT,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(h_LZvsSpectdiff_DMT,0)
%% entropy and Lz complexity
for s = 1:14
    [h_spect p_spect] = corrcoef(regDMT_min(:,s)',entsum_norm_dmt_min(:,s)');
    h_LZvsEnt_DMT(s) = h_spect(1,2);
    p_LZvsEnt_DMT(s) = p_spect(1,2);

end
%% plotting

figure
subplot(4,1,[1 2 3])

barh(1:14,h_LZvsEnt_DMT,'k')
hold on

set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if p_LZvsEnt_DMT(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(h_LZvsEnt_DMT(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(h_LZvsEnt_DMT(1,net)),' Pval: ',num2str(p_LZvsEnt_DMT(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(h_LZvsEnt_DMT,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(h_LZvsEnt_DMT,0)
%% correlations

figure
for s=1:14
subplot(3,5,s)
x3 = regDMT_min(:,s); % different data but the samee mean5HT2A_aalsymm(:,1); % taking the symmetric version
y3 = entsum_norm_dmt_min(:,s);

% correlation
[RHO,PVAL] = corr(x3,y3,'Type','Spearman')

p = polyfit(x3,y3,1);
% Evaluate the fitted polynomial p and plot:
f = polyval(p,x3);
plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'Intensity Ratings'})
ylabel('Entropy')
title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
set(gca, 'FontSize',14)

end

%% correlations

figure
for s=1:14
subplot(3,5,s)
x3 = regDMT_min(:,s); % different data but the same mean5HT2A_aalsymm(:,1); % taking the symmetric version
y3 = meas_min(:,s);

% correlation
[RHO,PVAL] = corr(x3,y3,'Type','Spearman')

p = polyfit(x3,y3,1);
% Evaluate the fitted polynomial p and plot:
f = polyval(p,x3);
plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
xlabel({'Intensity Ratings'})
ylabel('Entropy')
title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
set(gca, 'FontSize',14)

end
%% correlations

figure
for s=1:14
    subplot(3,5,s)
    x3 = regDMT(s,:)'; % different data but the same mean5HT2A_aalsymm(:,1); % taking the symmetric version
    y3 = meas(s,:)';
    
    % correlation
    [RHO,PVAL] = corr(x3,y3,'Type','Spearman')
    
    p = polyfit(x3,y3,1);
    % Evaluate the fitted polynomial p and plot:
    f = polyval(p,x3);
    plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
    xlabel({'Intensity Ratings'})
    ylabel('Entropy')
    title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
    set(gca, 'FontSize',14)

end
%% correlations

figure
for s=1:14
    subplot(3,5,s)
    x3 = RegDMT2(s,:)'; % different data but the same mean5HT2A_aalsymm(:,1); % taking the symmetric version
    y3 = squeeze(entsum_norm(1,s,:));
    
    % correlation
    [RHO,PVAL] = corr(x3,y3,'Type','Spearman')
    
    p = polyfit(x3,y3,1);
    % Evaluate the fitted polynomial p and plot:
    f = polyval(p,x3);
    plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
    xlabel({'Intensity Ratings'})
    ylabel('Entropy')
    title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
    set(gca, 'FontSize',14)

end