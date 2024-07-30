%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2004 Selen Atasoy
% Author: Jakub Vohryzek (jakub.vohryzek@upf.edu)
% Original script - author: Selen Atasoy (selenatasoy@gmail.com)

% Date: 20/11/2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_now   = 0;

%% folders & files
inputFile       = '/Users/jakub/Datasets/Psychedelics/DMT_cont/CH_Projections';
inputFile_CH    = '/Users/jakub/Datasets/Connectomes';
file_CH         = 'CH_sym_average_weighted.mat'; % original
file_dmt        = 'DMT_Projections_sym_prob_DOT.mat'; % original

%% load the files  
% load connectome harmonics
load(fullfile(inputFile_CH, file_CH));
%CH =  Connectome_Harmonics;
%V  = Eigenvalues(1:10000,1); % taking the first 100000 values
%% load dmt data
dmt = load('/Users/jakub/Datasets/Psychedelics/DMT_cont/CH_Projections/DMT_Projections_sym_prob_DOT.mat', 'Projections')
dmt = load(fullfile(inputFile, file_dmt));
mainNames_dmt = fieldnames(dmt.Projections);

%% define the dimensions
numSbj = size(dmt.Projections.(mainNames_dmt{1}), 1); % dmt number of subjects
numTp  = size(dmt.Projections.(mainNames_dmt{1}), 2); % dmt number of timepoints
numCH  = size(dmt.Projections.(mainNames_dmt{1}), 3); % number of harmonics
for f = 1:length(mainNames_dmt)
    dmt.Projections.(mainNames_dmt{f})     = dmt.Projections.(mainNames_dmt{f});
end
%% compute the energy and power distributions
% calculating Sdist for 15 bins across CH space
numBins = 15;
[CHBins,CHBinCenters] = surfFMRI_binCHinds(CH(1,1:numCH), numBins, 1); % 1 for logspace, 0 for linspace
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

entStats = squeeze(mean(entsum_norm,3));
%%
% figure
% yyaxis left
% plot(30:30:840,mean(intensityRatings_dmt_cont,1)','LineWidth',2)
% hold on
% plot(30:30:840,mean(intensityRatings_pcb_cont,1)','LineWidth',2)
% ylabel('Intensity Rating')
% yyaxis right
% plot(squeeze(nanmean(entsum_norm(1,:,:),2))','m-','LineWidth',2)
% plot(squeeze(nanmean(entsum_norm(2,:,:),2))','g-','LineWidth',2)
% plot([240 240],[0.82 0.91],'r-','LineWidth',2)
% xlim([0 840])
% set(gca,"XTick",0:30:840,'XTickLabel',0:28,'ycolor','m')
% legend('DMT intensity rating','PCB intensity rating','DMT','PCB','injection')
% ylabel('Entropy');xlabel('Minutes')

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
%ylim([0.84 0.89])
plot([8 8],[0.84 0.89],'r-','LineWidth',2)
xlim([0 28]);ylim([0.84 0.89])
set(gca,"XTick",0:1:28,'XTickLabel',0:1:28,'ycolor','m')
legend('DMT intensity rating','PCB intensity rating','DMT','PCB','injection')
ylabel('Entropy');xlabel('Minutes')
%%
for s = 1:14
    [h_spect, p_spect] = corr(intensityRatings_dmt_cont(s,:)',entsum_norm_dmt_min(:,s),'Type','Spearman')
    %tmp_PCB(s) = corrcoef(intensityRatings_pcb_cont(s,:)+0.01,entsum_norm_pcb_min(:,s)')
    h_meas_DMT(s) = h_spect;
    p_meas_DMT(s) = p_spect;

end

%% plotting
Sbj_names={'Sbj1','Sbj2','Sbj3','Sbj4','Sbj5','Sbj6','Sbj7',...
    'Sbj8','Sbj9','Sbj10','Sbj11','Sbj12','Sbj13','Sbj14'};

figure
subplot(4,1,[1 2 3])

barh(1:14,h_meas_DMT,'k')%,'FaceColor',color_vecs(net,:))
hold on

set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
title('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
xlim([0 1])

box off
for net=1:14
    if p_meas_DMT(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(h_meas_DMT(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(h_meas_DMT(1,net)),' Pval: ',num2str(p_meas_DMT(1,net)))};
    end
end

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(h_meas_DMT,'Group Average')
camroll(270)
ylim([0 1])
set(gca,'YTickLabel',[]);

[h p] = ttest(h_meas_DMT,0)
%% correlations
% 
% figure
% for s=1:14
% subplot(3,5,s)
% x3 = intensityRatings_dmt_cont(s,:)'; % different data but the samee mean5HT2A_aalsymm(:,1); % taking the symmetric version
% y3 = entsum_norm_dmt_min(:,s);
% 
% % correlation
% [RHO,PVAL] = corr(x3,y3,'Type','Spearman')
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
%%
% figure
% subplot(121)
% imagesc(tmp_DMT,[0 0.2]);axis square;colorbar
% subplot(122)
% imagesc(tmp_PCB,[0 0.2]);axis square;colorbar
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
figure,
plot(regDMT_min,'k');hold on;
plot(mean(regDMT_min,2),'y','LineWidth',2);hold on;
plot(RegDMT2_min,'k');hold on;
plot(mean(RegDMT2_min,2),'r','LineWidth',2);hold on;

%% 
figure
yyaxis left
plot(1:1:28,mean(Regdiff_min,2)','k-','LineWidth',2)
hold on
plot(1:1:28,mean(RegDMT2_min,2)','k-','LineWidth',2)
ylabel('LZ complexity')
yyaxis right
plot(squeeze(nanmean(entsum_norm_dmt_min,2))','m-','LineWidth',2)
plot(squeeze(nanmean(entsum_norm_pcb_min,2))','m--','LineWidth',2)
%ylim([0.84 0.89])
plot([8 8],[0.84 0.89],'r-','LineWidth',2)
xlim([0 28]);ylim([0.84 0.89])
set(gca,"XTick",0:1:28,'XTickLabel',0:1:28,'ycolor','m')
legend('Regdiff','RegDMT2','DMT CH repertoir entropy','PCB CH repertoir entropy','intervention')
ylabel('Entropy');xlabel('Minutes')
%%
[comp_entvsLZ(1,1) p_comp_entvsLZ(1,1)] = corr(mean(Regdiff_min,2),squeeze(nanmean(entsum_norm_dmt_min,2)),'type','Spearman')
[comp_entvsLZ(1,2) p_comp_entvsLZ(1,2)] = corr(mean(Regdiff_min,2),squeeze(nanmean(entsum_norm_pcb_min,2)),'type','Spearman')
[comp_entvsLZ(2,1) p_comp_entvsLZ(2,1)] = corr(mean(RegDMT2_min,2),squeeze(nanmean(entsum_norm_dmt_min,2)),'type','Spearman')
[comp_entvsLZ(2,2) p_comp_entvsLZ(2,2)] = corr(mean(RegDMT2_min,2),squeeze(nanmean(entsum_norm_pcb_min,2)),'type','Spearman')
% figure
% imagesc(comp_entvsLZ);colorbar
% set(gca,'XTickLabel',{'x','y'})

figure
imagesc(comp_entvsLZ)
for j1 = 1:2
    for j2 = 1:2
        %if j1 ~= j2
            caption = sprintf('%.2f', comp_entvsLZ(j1,j2));
            text(j2, j1, caption, 'FontSize', 15, 'Color', [0 0 0],'FontWeight','bold');
        %end
   end
end
colorbar

set(gca,'FontSize',14,'FontWeight','bold')
colormap('reds')
yticks(1:1:2); yticklabels({'Regdiff','RegDMT2'})
xticks(1:1:2); xticklabels({'DMT CH repertoir entropy','PCB CH repertoir entropy'})
xtickangle(45);ytickangle(45)
axis square

%% entropy and LZ complexity
for s = 1:14
    [h_LZvsEnt_DMT(s) p_LZvsEnt_DMT(s)] = corr(regDMT_min(:,s),entsum_norm_dmt_min(:,s),'Type','Spearman')
    [h_LZvsEnt_DMT2(s) p_LZvsEnt_DMT2(s)] = corr(RegDMT2_min(:,s),entsum_norm_dmt_min(:,s),'Type','Spearman')
    [h_LZvsEnt_DMTregdif(s) p_LZvsEnt_DMTregdif(s)] = corr(regdiff_min(:,s),entsum_norm_dmt_min(:,s),'Type','Spearman')
    [h_LZvsEnt_DMTRegdif(s) p_LZvsEnt_DMTRegdif(s)] = corr(Regdiff_min(:,s),entsum_norm_dmt_min(:,s),'Type','Spearman')

end
%% plotting - here change for the different regressors accordingly

figure
subplot(4,1,[1 2 3])
tmp_h = h_LZvsEnt_DMT;
tmp_p = p_LZvsEnt_DMT;

barh(1:14,tmp_h,'k')
hold on
set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if tmp_p(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(tmp_h(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(tmp_h(1,net)),' Pval: ',num2str(tmp_p(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end
title('LZvsEnt_DMT','Interpreter','none')

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(tmp_h,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(tmp_h,0)
%% plotting - here change for the different regressors accordingly

figure
subplot(4,1,[1 2 3])
tmp_h = h_LZvsEnt_DMT2;
tmp_p = p_LZvsEnt_DMT2;

barh(1:14,tmp_h,'k')
hold on
set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if tmp_p(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(tmp_h(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(tmp_h(1,net)),' Pval: ',num2str(tmp_p(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end
title('LZvsEnt_DMT2','Interpreter','none')

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(tmp_h,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(tmp_h,0)
%% plotting - here change for the different regressors accordingly

figure
subplot(4,1,[1 2 3])
tmp_h = h_LZvsEnt_DMTregdif;
tmp_p = p_LZvsEnt_DMTregdif;

barh(1:14,tmp_h,'k')
hold on
set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if tmp_p(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(tmp_h(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(tmp_h(1,net)),' Pval: ',num2str(tmp_p(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end
title('LZvsEnt_DMTregdif','Interpreter','none')

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(tmp_h,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(tmp_h,0)
%% plotting - here change for the different regressors accordingly

figure
subplot(4,1,[1 2 3])
tmp_h = h_LZvsEnt_DMTRegdif;
tmp_p = p_LZvsEnt_DMTRegdif;
barh(1:14,tmp_h,'k')
hold on
set(gca,'YTick',1:14);set(gca,'YTickLabel',Sbj_names,'Fontsize',8)
set(gca,'YTickLabelRotation',45)
xlabel('Correlation','Fontsize',16)% ;ylabel('Correlation','Fontsize',16)
%ylim([-.4 1.1])
box off
for net=1:14
    if tmp_p(1,net)<(0.05)
        % To plot an asterisk when significant
        plot(tmp_h(1,net)*1.1,net,'*k')
            string2annot = {strcat('Corr: ',num2str(tmp_h(1,net)),' Pval: ',num2str(tmp_p(1,net)))};
        %text(5,0.4,string2annot,'FontSize',20)
    end
end
title('LZvsEnt_DMTRegdif','Interpreter','none')

set(gca,'YGrid','on','FontSize',14)
subplot(4,1,4)

boxplot(tmp_h,'Group Average')
camroll(270)
%ylim([-0.3 0.7])
set(gca,'YTickLabel',[]);

[h p] = ttest(tmp_h,0)
%%
% %% correlations
% 
% figure
% for s=1:14
% subplot(3,5,s)
% x3 = regDMT_min(:,s); % different data but the samee mean5HT2A_aalsymm(:,1); % taking the symmetric version
% y3 = entsum_norm_dmt_min(:,s);
% 
% % correlation
% [RHO,PVAL] = corr(x3,y3,'Type','Pearson')
% 
% p = polyfit(x3,y3,1);
% % Evaluate the fitted polynomial p and plot:
% f = polyval(p,x3);
% plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
% xlabel({'LZc'})
% ylabel('Entropy')
% title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
% set(gca, 'FontSize',14)
% 
% end
% %% correlations
% 
% figure
% for s=1:14
%     subplot(3,5,s)
%     x3 = RegDMT2(s,:)'; % different data but the same mean5HT2A_aalsymm(:,1); % taking the symmetric version
%     y3 = squeeze(entsum_norm(1,s,:));
% 
%     % correlation
%     [RHO,PVAL] = corr(x3,y3,'Type','Spearman')
% 
%     p = polyfit(x3,y3,1);
%     % Evaluate the fitted polynomial p and plot:
%     f = polyval(p,x3);
%     plot(x3,y3,'k.',x3,f,'r-','MarkerSize',20,'LineWidth',4);axis square
%     xlabel({'LZc'})
%     ylabel('Entropy')
%     title({['Sub #',num2str(s),' Spearman corr: ' num2str(round(RHO,3))],['p-value of ' num2str(round(PVAL,3))]})
%     set(gca, 'FontSize',14)
% 
% end