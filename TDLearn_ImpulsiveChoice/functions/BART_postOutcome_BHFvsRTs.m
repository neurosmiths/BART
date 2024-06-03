function [] = BART_postOutcome_BHFvsRTs(ptID)
% BART_BALOONONSET_POSTOUTCOME_BHF examine post outcome encoding of risk/reward
%
%   [BARTstats] = BART_balloonOnset_postOutcome_bhf(ptID,nevFile) analyzes LFP data
%   for the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20200823

ptID = '201901'; % this line only for debugging


BARTdir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data'];

nevList = dir([BARTdir '\*.nev']);

% loading PTB BHV
load(fullfile(BARTdir,[ptID '.bartBHV.mat']))

% electrode information.
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

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

% timing parameters.
pre = 2;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% % task parameters (out of order)
% alignName = 'banksandpops';
% bankTimes = trigTimes(trigs==25);
% popTimes = trigTimes(trigs==26);
% nBanks = length(bankTimes);
% nPops = length(popTimes);
% outcomeType = [ones(1,nBanks) 2*ones(1,nPops)]; % 1 = bank, 2 = pop
% nTrials = length(outcomeType);
% %1: trial start ::      [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]

% task parameters in chronological order..
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% only including complete trials.
nTrials = length(respTimes);
balloonIDs = balloonIDs(1:nTrials);
balloonType = balloonType(1:nTrials);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
isCTRL = isCTRL(1:nTrials);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
postOutcome = [NaN; outcomeType(1:end-1)]; %

% colormap for balloon color. #uselelss???
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% outcome color map
oMap(1,:) = rgb('forestgreen');
oMap(2,:) = rgb('orangered');

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% high notch for one patient with a high noise harmonic
if isequal(ptID,'201903')
    [b120,a120] = iirnotch(120/(Fs/2),60/(Fs/2)/35);
end

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% theta filter
[bTh,aTh] = butter(4,[4 8]/(Fs/2));

% smoothing
smoothType = 'movmean';
smoothFactor = Fs./5;

% aligning on balloon appearances to start.
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    CUEmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    OUTmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        % matrix of high gamma signals. [channels X samples X trtials]
        cueHG(ch,:,tt) = smoothdata(abs(hilbert(filtfilt(b,a,double(CUEmat(ch,:,tt)))))',smoothType,smoothFactor);
        outHG(ch,:,tt) = smoothdata(abs(hilbert(filtfilt(b,a,double(OUTmat(ch,:,tt)))))',smoothType,smoothFactor);

        % theta power
        cueTHETA(ch,:,tt) = smoothdata(abs(hilbert(filtfilt(bTh,aTh,double(CUEmat(ch,:,tt)))))',smoothType,smoothFactor);
        outTHETA(ch,:,tt) = smoothdata(abs(hilbert(filtfilt(bTh,aTh,double(OUTmat(ch,:,tt)))))',smoothType,smoothFactor);
    end
end


%% reaction time stuff
RTs = [data(1:nTrials).rt]';
WOI = (tSec>0.25 & tSec<1.25);
outlyingRTs = outliers(RTs);

%% points
pointsEarned = [data(1:nTrials).points];


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
expectedOutcome(balloonIDs==1) = cumsum(outcomeType(balloonIDs==1)==1)./(1:sum(balloonIDs==1))';
expectedOutcome(balloonIDs==2) = cumsum(outcomeType(balloonIDs==2)==1)./(1:sum(balloonIDs==2))';
expectedOutcome(balloonIDs==3) = cumsum(outcomeType(balloonIDs==3)==1)./(1:sum(balloonIDs==3))';

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

% previous trial indices
prevTrialRiskGroups = [-10; riskGroups(1:end-1)];
prevTrialOutcome = [-10; outcomeType(1:end-1)];
prevTrialCTRL = [-10; isCTRL(1:end-1)];
prevTrialPop_currentTrialBank = prevTrialOutcome==2 & outcomeType==1;
prevTrialBank_currentTrialPop = prevTrialOutcome==1 & outcomeType==2;
prevTrialPop_currentTrialPop = prevTrialOutcome==2 & outcomeType==2;
prevTrialBank_currentTrialBank = prevTrialOutcome==1 & outcomeType==1;

% behavioral figure
figure
hold on
violinPlot(RTs(~prevTrialCTRL & prevTrialPop_currentTrialPop),'showMM',4,'color',rgb('lime'),'xValues',1)
violinPlot(RTs(~prevTrialCTRL & prevTrialBank_currentTrialPop),'showMM',4,'color',rgb('aqua'),'xValues',2)
violinPlot(RTs(~prevTrialCTRL & prevTrialPop_currentTrialBank),'showMM',4,'color',rgb('fuchsia'),'xValues',3)
violinPlot(RTs(~prevTrialCTRL & prevTrialBank_currentTrialBank),'showMM',4,'color',rgb('red'),'xValues',4)
hold off


% looping over channels
for ch2 = 1:nChans

        % plotting
        if ishandle(ch2); close(ch2); end
        figure(ch2)

        %% plotting the signals.
        % cue
        subplot(2,3,1)
        hold on
        patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialPop),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialPop),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialPop_currentTrialPop)))))...
            fliplr(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialPop),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialPop),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialPop_currentTrialPop)))))]...
            ,rgb('lime'),'facealpha',0.3,'edgecolor','none')
        patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialPop),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialPop),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialBank_currentTrialPop)))))...
            fliplr(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialPop),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialPop),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialBank_currentTrialPop)))))]...
            ,rgb('aqua'),'facealpha',0.3,'edgecolor','none')
        patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialBank),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialBank),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialPop_currentTrialBank)))))...
            fliplr(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialBank),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialBank),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialPop_currentTrialBank)))))]...
            ,rgb('fuchsia'),'facealpha',0.3,'edgecolor','none')
        patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialBank),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialBank),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialBank_currentTrialBank)))))...
            fliplr(mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialBank),3)-(std(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialBank),[],3)./(sqrt(sum(~prevTrialCTRL & prevTrialBank_currentTrialBank)))))]...
            ,rgb('red'),'facealpha',0.3,'edgecolor','none')
        plot(tSec,mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialPop),3),'color',rgb('lime'))
        plot(tSec,mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialPop),3),'color',rgb('aqua'))
        plot(tSec,mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialPop_currentTrialBank),3),'color',rgb('fuchsia'))
        plot(tSec,mean(cueHG(ch2,:,~prevTrialCTRL & prevTrialBank_currentTrialBank),3),'color',rgb('red'))

        % deets
        hold off
        axis tight square
        xlim([-pre+1 post-1])
        xlabel('time relative to balloon onset (s)')
        ylabel('broadband high frequency power')


        % saving figures
        saveDir = 'D:\Data\Rhiannon\BART_postOutcome\'; %does not exist
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        halfMaximize(ch2,'left')
        saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_postOutcomeXcurrentOutcome_balloonOnset_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))

end









 
% 
% %%
% keyboard
% 
% 
%     %     cueLM = fitlm(RTs,squeeze(mean(cueHG(ch2,WOI,:),2)),'Exclude',outlyingRTs);
%     %     outLM = fitlm(RTs,squeeze(mean(outHG(ch2,WOI,:),2)),'Exclude',outlyingRTs);
%     %     postOutLM = fitlm(RTs(1:end-1),squeeze(mean(cueHG(ch2,WOI,2:end),2)),'Exclude',outlyingRTs(1:end-1));
% 
%     crit = 0.05;
%     if (cueLM.Coefficients.pValue(2)<crit || outLM.Coefficients.pValue(2)<crit || postOutLM.Coefficients.pValue(2)<crit)
% 
% 
%     end
% 
% 
% hold on
% for rsk = 1:4
%     patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,riskGroups==(rsk-1)),3)-(std(cueHG(ch2,:,riskGroups==(rsk-1)),[],3)./(sqrt(sum(riskGroups==(rsk-1))))))...
%         fliplr(mean(cueHG(ch2,:,riskGroups==(rsk-1)),3)-(std(cueHG(ch2,:,riskGroups==(rsk-1)),[],3)./(sqrt(sum(riskGroups==(rsk-1))))))]...
%         ,cMap(rsk,:),'facealpha',0.3,'edgecolor','none')
%     plot(tSec,mean(cueHG(ch2,:,riskGroups==(rsk-1)),3),'color',cMap(rsk,:))
% end
% hold off
% 
% % deets
% hold off
% axis tight square
% xlim([-pre+1 post-1])
% xlabel('time relative to balloon onset (s)')
% ylabel('broadband high frequency power')
% 
% 
% % plotting the risk and reward categories
% subplot(2,3,2)
% hold on
% for rsk = 1:7
%     patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,rewardGroups==(rsk-1)),3)-(std(cueHG(ch2,:,rewardGroups==(rsk-1)),[],3)./nTrials))...
%         fliplr(mean(cueHG(ch2,:,rewardGroups==(rsk-1)),3)-(std(cueHG(ch2,:,rewardGroups==(rsk-1)),[],3)./nTrials))]...
%         ,rcMap(rsk,:),'facealpha',0.3,'edgecolor','none')
%     plot(tSec,mean(cueHG(ch2,:,rewardGroups==(rsk-1)),3),'color',rcMap(rsk,:))
% end
% hold off
% 
% % deets
% hold off
% axis tight square
% xlim([-pre+1 post-1])
% xlabel('time relative to balloon onset (s)')
% ylabel('broadband high frequency power')
% try
%     title(anatomicalLocs(ch2))
% catch
%     title('uh oh... no anatomical label...')
% end
% 
% subplot(2,3,3)
% hold on
% for rsk = 1:4
%     patch([tSec fliplr(tSec)],[(mean(cueHG(ch2,:,prevTrialRiskGroups==(rsk-1)),3)-(std(cueHG(ch2,:,prevTrialRiskGroups==(rsk-1)),[],3)./(sqrt(sum(prevTrialRiskGroups==(rsk-1))))))...
%         fliplr(mean(cueHG(ch2,:,prevTrialRiskGroups==(rsk-1)),3)-(std(cueHG(ch2,:,prevTrialRiskGroups==(rsk-1)),[],3)./(sqrt(sum(prevTrialRiskGroups==(rsk-1))))))]...
%         ,cMap(rsk,:),'facealpha',0.3,'edgecolor','none')
%     plot(tSec,mean(cueHG(ch2,:,prevTrialRiskGroups==(rsk-1)),3),'color',cMap(rsk,:))
% end
% hold off
% 
% % deets
% hold off
% axis tight square
% xlim([-pre+1 post-1])
% xlabel('time relative to balloon onset (s)')
% ylabel('broadband high frequency power')
% title('colored by previous trial cue')
% 
% 
% %% plotting the models
% % cue-RT
% subplot(2,3,4)
% h = plotAdded(cueLM);
% h(1).Marker = '.';
% h(1).Color = 'k';
% legend off
% xlabel('RT')
% ylabel('cue-aligned HFA')
% axis tight square
% title(sprintf('t(%d) = %.2f, p = %.2f',nTrials-2,cueLM.Coefficients.tStat(2),cueLM.Coefficients.pValue(2)))
% 
% 
% % outcome-RT
% subplot(2,3,5)
% h = plotAdded(outLM);
% h(1).Marker = '.';
% h(1).Color = 'k';
% legend off
% xlabel('RT')
% ylabel('outcome-aligned HFA')
% axis tight square
% title(sprintf('t(%d) = %.2f, p = %.2f',nTrials-2,outLM.Coefficients.tStat(2),outLM.Coefficients.pValue(2)))
% 
% % post-outcome-cue RT
% subplot(2,3,6)
% h = plotAdded(postOutLM);
% h(1).Marker = '.';
% h(1).Color = 'k';
% legend off
% xlabel('RT')
% ylabel('post-outcome cue-aligned HFA')
% axis tight square
% title(sprintf('t(%d) = %.2f, p = %.2f',nTrials-2,postOutLM.Coefficients.tStat(2),postOutLM.Coefficients.pValue(2)))

