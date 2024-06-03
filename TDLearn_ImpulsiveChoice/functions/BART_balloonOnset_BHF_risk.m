function [BARTstats] = BART_balloonOnset_BHF_risk(ptID,BARTdir)
% BART_baloonOnset_BHF_risk analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_balloonOnset_BHF_risk(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   [BARTstats] = BART_balloonOnset_BHF_risk(ptID,nevFile,chan2Plot) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.

% author: EHS20181005

% ptID = '202009'; % this line for debugging only

% getting data dir
if nargin<2
    BARTdir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
end
nevList = dir([BARTdir '/*.nev']);

% getting complete behavioral data
bhvList = dir([BARTdir '/*BHV.mat'])
load(fullfile(bhvList.folder,bhvList.name))

% getting electrode information
[trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);

% load and define triggers from nevFle
if length(nevList)>2
    error('2manyNEVs!')
else
    nevFile = fullfile(BARTdir,nevList(1).name);
end
NEV = openNEV(fullfile(BARTdir,nevList(1).name),'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
nChans = length(trodeLabels(isECoG));
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% task timings
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;
pointsEarned = [data.points];

% only including complete trials; generally => excluding the last trial.
nTrials = min([length(respTimes) length(balloonIDs)]);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
isCTRL = isCTRL(1:nTrials);
pointsEarned = pointsEarned(1:nTrials);
poppedTrials = logical(trigs(trigs==25 | trigs==26)-25); % 0 = bank; 1 = pop;


% duration of balloon inflation
inflateDurations = respTimes-inflateTimes;

% timing parameters.
pre = 2;
post = 4;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% risk colormap
cMap(2,:) = [1 0.9 0];
cMap(3,:) = [1 0.5 0];
cMap(4,:) = [1 0 0];
cMap(1,:) = [0.5 0.5 0.5];

% risk and reward colormap
rcMap(1,:) = [0.5 0.5 0.5];
rcMap(2,:) = [1 0 0];
rcMap(3,:) = [1 0.5 0];
rcMap(4,:) = [1 0.9 0];
rcMap(5,:) = rgb('lightcoral');
rcMap(6,:) = rgb('rosybrown');
rcMap(7,:) = rgb('violet');


% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% aligning on balloon appearances to start.
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        % matrix of high gamma signals. [channels X samples X trtials]
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
end


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


% how much do we smooth the data?
% Fs = 1000 for these experiments.
% so this would smooth by 200 ms moving average.
smoothFactor = Fs./5;
% plotting responses
for ch2 = 1:nChans
    % setting up variables for statistics
    statTimeWin = [0 1]; % in seconds.
    riskVar = squeeze(mean(HGmat(ch2,tSec>statTimeWin(1) & tSec<statTimeWin(2),:),2));
    
    %% anovas for processing of risk, while controlling for reward.  incrasing risk => gray, yellow, orange, red
    % risk
    [BARTstats(ch2).pRisk,BARTstats(ch2).tblRisk,BARTstats(ch2).statsRisk] = anova1(riskVar,riskGroups);
    
    % risk controlling for reward - 2 way anova
    % TODO:: this model is singular... need a better way of doing this.
    % maybe with a continuous baloon size variable?
    rskTbl = table(riskVar,riskGroups,rewardGroups,expectedReward,expectedOutcome,'VariableNames',{'HG','riskCats','rewardCats','expectedReward','expectedOutcome'});
    % glme = fitglme(rskTbl,'HG ~ 1 + expectedReward*riskCats');
    formula = 'HG ~ 1 + expectedReward*expectedOutcome + (1 + HG | riskCats) + (1 + HG | rewardCats)';
    glme = fitglme(rskTbl,formula);
    BARTstats(ch2).statsRiskReward = glme.anova;
    close all
    
    % [20200812] There was a bug in my definition of rewardVar so this may
    % work now...
    % BARTstats(ch2).pRiskReward,BARTstats(ch2).tblRiskReward,BARTstats(ch2).statsRiskReward] = anovan(riskVar,{riskGroups,rewardGroups},'model','interaction','varnames',{'potentialRisk','potentialReward'});
    
    % plotting
    if ishandle(ch2); close(ch2); end
    figure(ch2)
    
    subplot(3,2,1)
    hold on
    for rsk = 1:4
        patch([tSec fliplr(tSec)],[(smooth(mean(HGmat(ch2,:,riskGroups==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,riskGroups==(rsk-1)),[],3),smoothFactor)'./nTrials))...
            fliplr(smooth(mean(HGmat(ch2,:,riskGroups==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,riskGroups==(rsk-1)),[],3),smoothFactor)'./nTrials))]...
            ,cMap(rsk,:),'facealpha',0.3,'edgecolor','none')
        plot(tSec,smooth(mean(HGmat(ch2,:,riskGroups==(rsk-1)),3),smoothFactor)','color',cMap(rsk,:))
    end
    hold off
    
    % deets
    hold off
    axis tight square
    xlim([-pre+1 post-1])
    xlabel('time relative to balloon onset (s)')
    ylabel('broadband high frequency power')
    title(sprintf('%s (%s) 1-way ANOVA: p = %.2f, t(%d) = %.2f',...
        deblank(trodeLabels{ch2}),char(anatomicalLocs{ch2}),BARTstats(ch2).pRisk,BARTstats(ch2).tblRisk{2,3},BARTstats(ch2).tblRisk{2,5}))
    
    % plotting the risk and reward categories
    subplot(3,2,2)
    hold on
    for rsk = 1:7
        patch([tSec fliplr(tSec)],[(smooth(mean(HGmat(ch2,:,rewardGroups==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,rewardGroups==(rsk-1)),[],3),smoothFactor)'./nTrials))...
            fliplr(smooth(mean(HGmat(ch2,:,rewardGroups==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,rewardGroups==(rsk-1)),[],3),smoothFactor)'./nTrials))]...
            ,rcMap(rsk,:),'facealpha',0.3,'edgecolor','none')
        plot(tSec,smooth(mean(HGmat(ch2,:,rewardGroups==(rsk-1)),3),smoothFactor)','color',rcMap(rsk,:))
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
    
    % saving statistics.
    save(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/%s_riskAndRewardANOVAs.mat',ptID,ptID),'BARTstats','-v7.3')
    
    
end % looping over channels

end % eof





