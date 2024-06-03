function [] = BART_balloonOnset_FR_risk(ptID)   % BARTstats
% BART_BALLOONONSET_FR_RISK analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_balloonOnset_FR_risk(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005

% input args
ptID = '202107';

% loading data.
nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
% [trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% loading behavioral matFile
matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)
pointsEarned = [data.points];

% task parameters in chronological order..
% standard 3-D [chan, unit, timestamp (seconds)] matrix.
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets.
inclChans = unique(ChanUnitTimestamp(:,1));
microLabels = microLabelsBART(ptID);
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.\
nChans = length(inclChans);

% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% adjusting for trial numbers
% only including complete trials; generally => excluding the last trial.
nTrials = min([length(respTimes) length(balloonIDs)]);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
isCTRL = isCTRL(1:nTrials);
pointsEarned = pointsEarned(1:nTrials);
poppedTrials = logical(trigs(trigs==25 | trigs==26)-25); % 0 = bank; 1 = pop;
[~,sortedTrialIdcs] = sort(balloonIDs);

% duration of balloon inflation
inflateDurations = respTimes-inflateTimes;

% timing parameters.
pre = 3;
post = 3;


%% color maps
% risk colormap
cMap(2,:) = [1 0.9 0];
cMap(3,:) = [1 0.5 0];
cMap(4,:) = [1 0 0];
cMap(1,:) = [0.5 0.5 0.5];

% risk and reward colormap
rcMap(1,:) = [0.5 0.5 0.5];
rcMap(2,:) = [1 0 0];
rcMap(3,:) = [1 0.5 0];
racMap(4,:) = [1 0.9 0];
rcMap(5,:) = rgb('lightcoral');
rcMap(6,:) = rgb('rosybrown');
rcMap(7,:) = rgb('violet');


%% setting up regressors for linear models.
% [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]

% making a variable that represents cumulative
expectedReward = zeros(nTrials,1);
expectedReward(balloonIDs==1 | balloonIDs==11) = cumsum(pointsEarned(balloonIDs==1 | balloonIDs==11))./(1:sum(balloonIDs==1 | balloonIDs==11));
expectedReward(balloonIDs==2 | balloonIDs==12) = cumsum(pointsEarned(balloonIDs==2 | balloonIDs==12))./(1:sum(balloonIDs==2 | balloonIDs==12));
expectedReward(balloonIDs==3 | balloonIDs==13) = cumsum(pointsEarned(balloonIDs==3 | balloonIDs==13))./(1:sum(balloonIDs==3 | balloonIDs==13));
expectedReward(balloonIDs==14) = 0;

% making a variable containing expected outcome
expectedOutcome = double(isCTRL);
expectedOutcome(balloonIDs==1) = cumsum(~poppedTrials(balloonIDs==1))./(1:sum(balloonIDs==1))';
expectedOutcome(balloonIDs==2) = cumsum(~poppedTrials(balloonIDs==2))./(1:sum(balloonIDs==2))';
expectedOutcome(balloonIDs==3) = cumsum(~poppedTrials(balloonIDs==3))./(1:sum(balloonIDs==3))';

% [20200807] risk variables - controls = zero. otherwise, identities remain the same.
riskGroups = balloonIDs;
riskGroups(isCTRL) = 0;

% [20200807] reward vraiables - gray = zero; red, orange, yellow; incremented red, orange and yellow visual controls, which may not be the best approach, but reasonable.
rewardGroups = balloonIDs;
rewardGroups(balloonIDs==14) = 0;
rewardGroups(balloonIDs==1) = 3;
rewardGroups(balloonIDs==3) = 1;
rewardGroups(balloonIDs==11) = 6;
rewardGroups(balloonIDs==12) = 5;
rewardGroups(balloonIDs==13) = 4;


% looping over Channels
for ch = 1:nChans
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));
    for un = 1:nUnits
        fprintf('\nprocessing and plotting for unit %d of %d',un,nUnits)
        
        % getting unit times for the current channel and unit.
        unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds

        % loooping over trials
        for tt = 1:nTrials
            
            % putting the data in a structure
            spikes.channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post) - repmat(balloonTimes(tt)-pre,length(unitTimes(unitTimes>balloonTimes(tt)-pre & unitTimes<balloonTimes(tt)+post)),1);
            
        end % looping over trials
        
        % calculating psths
        kernelWidth = 50  ./1000;
        [Rred,t,Ered] = psth(spikes.channel(ch).unit(un).trial(balloonIDs==1 | balloonIDs==11), kernelWidth, 'n', [0 pre+post]);
        [Rorange,t,Eorange] = psth(spikes.channel(ch).unit(un).trial(balloonIDs==2 | balloonIDs==12), kernelWidth, 'n', [0 pre+post]);
        [Ryellow,t,Eyellow] = psth(spikes.channel(ch).unit(un).trial(balloonIDs==3 | balloonIDs==13), kernelWidth, 'n', [0 pre+post]);
        [Rgray,t,Egray] = psth(spikes.channel(ch).unit(un).trial(balloonIDs==14), kernelWidth, 'n', [0 pre+post]);
        
        
        for t2 = 1:nTrials
            % single trial firing rates.
            if isempty(spikes.channel(ch).unit(un).trial(t2).times)
                Rst(t2,:) = zeros(1,length(t));
                fprintf('\n     zero spikes recorded in trial %d! \n', t2)
            else
                [Rst(t2,:)] = psth(spikes.channel(ch).unit(un).trial(t2), kernelWidth, 'n', [0 pre+post]);
            end
        end % looping over trials
        
        % timing
        tSec = t-repmat(pre,1,length(t));
        
        axis square
        if ch<=8
            title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{1}))
        elseif ch>8 & ch<=16
            title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{2}))
        elseif ch>16 & ch<=24
            title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{3}))
        elseif ch>24 & ch<=32
            title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,'4th BF???'))
        else
            error('trying to loop over more channels than recorded electrodes. something went wrong...')
        end
        xlim([-1 2])
        ylabel('trials')
        
        
        % setting up variables for statistics
        statTimeWin = [0 1]; % in seconds.
        riskVar = squeeze(mean(Rst(:,tSec>statTimeWin(1) & tSec<statTimeWin(2)),2));
        
        %% anovas for processing of risk, while controlling for reward.  incrasing risk => gray, yellow, orange, red
        % risk
        [BARTstats(ch,un).pRisk,BARTstats(ch,un).tblRisk,BARTstats(ch,un).statsRisk] = anova1(riskVar,riskGroups);
        
        % risk controlling for reward - 2 way anova
        % TODO:: this model is singular... need a better way of doing this.
        % maybe with a continuous baloon size variable?
        rskTbl = table(riskVar,riskGroups,rewardGroups,expectedReward,expectedOutcome,'VariableNames',{'HG','riskCats','rewardCats','expectedReward','expectedOutcome'});
        % glme = fitglme(rskTbl,'HG ~ 1 + expectedReward*riskCats');
        formula = 'HG ~ 1 + expectedReward*expectedOutcome + (1 + HG | riskCats) + (1 + HG | rewardCats)';
        glme = fitglme(rskTbl,formula);
        BARTstats(ch,un).statsRiskReward = glme.anova;
        close all
        
        % [20200812] There was a bug in my definition of rewardVar so this may
        % work now...
        % BARTstats(ch2).pRiskReward,BARTstats(ch2).tblRiskReward,BARTstats(ch2).statsRiskReward] = anovan(riskVar,{riskGroups,rewardGroups},'model','interaction','varnames',{'potentialRisk','potentialReward'});
        
        % plotting
        if ishandle(un); close(un); end
        figure(un)
        
        subplot(3,2,2)
        hold on
        for rsk = 1:4
            patch([tSec fliplr(tSec)],[((mean(Rst(riskGroups==(rsk-1),:)))'-((std(Rst(riskGroups==(rsk-1),:)))'./nTrials))'...
                fliplr((mean(Rst(riskGroups==(rsk-1),:)))'-((std(Rst(riskGroups==(rsk-1),:)))'./nTrials))']...
                ,cMap(rsk,:),'facealpha',0.3,'edgecolor','none')
            plot(tSec,(mean(Rst(riskGroups==(rsk-1),:)))','color',cMap(rsk,:))
        end
        hold off
        
        % deets
        hold off
        axis tight square
        xlim([-pre+1 post-1])
        xlabel('time relative to balloon onset (s)')
        ylabel('firing rates')
        title(sprintf('%s (%s) 1-way ANOVA: p = %.2f, t(%d) = %.2f',...
            deblank(trodeLabels{ch}),char(anatomicalLocs{ch2}),BARTstats(ch2).pRisk,BARTstats(ch2).tblRisk{2,3},BARTstats(ch2).tblRisk{2,5}))
        
        % plotting the risk and reward categories
        subplot(3,2,2)
        hold on
        for rsk = 1:7
            patch([tSec fliplr(tSec)],[((mean(Rst(rewardGroups==(rsk-1),:)))'-((std(Rst(rewardGroups==(rsk-1),:)))'./nTrials))'...
                fliplr((mean(Rst(rewardGroups==(rsk-1),:)))'-((std(Rst(rewardGroups==(rsk-1),:)))'./nTrials))']...
                ,rcMap(rsk,:),'facealpha',0.3,'edgecolor','none')
            plot(tSec,(mean(Rst(rewardGroups==(rsk-1),:)))','color',rcMap(rsk,:))
        end
        hold off
        
        % deets
        hold off
        axis tight square
        xlim([-pre+1 post-1])
        xlabel('time relative to balloon onset (s)')
        ylabel('broadband high frequency power')
        title(sprintf('%s (%s) LMM (%s) ANOVA interaction: p = %.2f, F(%d) = %.2f',...
            deblank(trodeLabels{ch2}),char(anatomicalLocs{ch2}),formula,BARTstats(ch2).statsRiskReward.pValue(4),BARTstats(ch2).statsRiskReward.DF2(4),BARTstats(ch2).statsRiskReward.FStat(4)))
        
        subplot(3,2,3)
        notBoxPlot(riskVar,riskGroups+1);
        % plotting pairwise comparisons for risk.
        if BARTstats(ch2).pRisk<0.05
            pairwiseComps = multcompare(BARTstats(ch2).statsRisk,'Display','off');
            sigComps = pairwiseComps(pairwiseComps(:,6)<0.05,:);
            hold on
            for ls = 1:sum(pairwiseComps(:,6)<0.05)
                line([sigComps(ls,1) sigComps(ls,1)],[max(riskVar)+0.4 max(riskVar)+0.4],'linewidth',2)
            end
        end
        axis square
        xlabel('increasing risk: none, yellow, orange, red')
        ylabel('mean BHF power')
        title('lines: significant pairwise comparisons')
        
        
        subplot(3,2,4)
        notBoxPlot(riskVar,rewardGroups+1)
        axis square
        xlabel('increasing potential reward: gray, active[yellow, orange, red], passive[yellow, orange, red]')
        ylabel('expectedReward')
        zlabel('linear model fixed effects (adjusted high gamma)')
        title('lines: significant pairwise comparisons')
        
        
        % saving figure
        saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s',ptID);
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        halfMaximize(ch2,'left')
        saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_balloonOnsetRisk_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        if BARTstats(ch2).pRisk<0.05
            saveas(ch2,fullfile('/home/user1/Dropbox/significantRiskTrodes','riskOneWay',sprintf('pt%s_%s_balloonOnsetRisk_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        end
        
        if BARTstats(ch2).statsRiskReward.pValue(4)<0.05
            saveas(ch2,fullfile('/home/user1/Dropbox/significantRiskTrodes','riskRewardTwoWay',sprintf('pt%s_%s_balloonOnsetRisk_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        end
        close all
        
        
        
        % saving figures
        saveDir = sprintf('~/Dropbox/BART_firing/');
        if ~isempty(Rred)
            %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
            saveas(un,fullfile(saveDir,sprintf('pt%s_ch%d_un%d_balloonOnset_firingrates.pdf',ptID,ch,un)))
        end
        close(un)
        
    end % looping over units
end % looping over channels



















%
%
%
% %% loading kilosorted spikes
% spikesDir = '/media/user1/data4TB/data/spikes_BART_202002';
% spikeStruct = loadKSdir(spikesDir);
% clusterIDs = unique(spikeStruct.clu);
% nUnits = length(clusterIDs);
%
%
% %% spike parameters
% binSize = 10./1000;
% kernelWidth = 50./1000;
%
%
% alignSpots = outcomeTimes;
% % epoch Firing Rates
% for un = 1:nUnits
%     thisClust = spikeStruct.clu;
%     unitTimes = (spikeStruct.st(spikeStruct.clu==clusterIDs(un))');
%
%     % compute everything
%     [psth, bins, rasterX, rasterY, ~, ba] = psthAndBA(unitTimes, alignSpots, [-pre post], binSize);
%
%     % PSTH smoothing filter
%     gw = gausswin(round(kernelWidth*600),3); % kernelwidth[*4]?? in units of ms here, because gausswin fuction likes ints.
%     smWin = gw./sum(gw);
%
%     % smooth ba
%     baSm = conv2(smWin,1,ba', 'same')'./binSize;
%
%     % stacking smoothed psths into a tensor
%     softNorm = false;
%     if softNorm % [20181114] this is what churchland does for all of his papes.
%         trialTensor(un,:,:) = baSm./(range(range(baSm))+5);
%     else
%         trialTensor(un,:,:) = baSm;
%     end
%     % [20181114] Better to observe raw firing rates to start.
%     baTensor(un,:,:) = ba;
%
%     % timing vector
%     tSec = linspace(-pre,post,size(trialTensor,3));
%
%     % visualize the data here for each unit
%     figure(un)
%     subplot(2,1,1)
%     plotSpikeRaster(logical(ba),'PlotType','vertline');
%     xlim([200 400])
%     axis square off
%     title(sprintf('unitID: %d',clusterIDs(un)))
%
%     subplot(2,1,2)
%     hold on
%     patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==2,:),2)-(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))'...
%         fliplr(squeeze(mean(trialTensor(un,outcomeType==2,:),2)+(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))')]...
%         ,rgb('orangered'),'facealpha',0.5,'edgecolor','none')
%     plot(tSec,squeeze(mean(trialTensor(un,outcomeType==2,:),2))','color',rgb('orangered'))
%     text(0,2.5,sprintf('N = %d',sum(outcomeType==2)),'color',rgb('orangered'))
%
%     patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==1,:),2)-(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))'...
%         fliplr(squeeze(mean(trialTensor(un,outcomeType==1,:),2)+(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))')]...
%         ,rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
%     plot(tSec,squeeze(mean(trialTensor(un,outcomeType==1,:),2))','color',rgb('forestgreen'))
%     text(0,3.5,sprintf('N = %d',sum(outcomeType==1)),'color',rgb('forestgreen'))
%     %     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
%     hold off
%     axis tight square
%     xlim([-1 2])
%     xlabel('time relative to outcome (s)')
%     ylabel('firing rate (spks/s)')
%     title(sprintf(''))
%     halfMaximize(un,'left')
%
%     % saving figures
%     saveDir = sprintf('~/Figs/BART/%s',ptID);
%     if exist(saveDir,'dir')
%         %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
%         saveas(un,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_firing rates.pdf',ptID,deblank(trodeLabels{ch2}))))
%     end
%
%
% end
%
%
%
