function [TDdata,rewardData,riskData,rewardAndRiskHGs] = BART_bankpop_bhf_TDlearn(ptID,whichRL,whichTrialsLM,whichTrialsTD,plotFlag)   % BARTstats
% BART_BANKPOP_BHF_TDLEARN fits BART responses to Q-learning model.
%
%   [TDdata] = BART_bankpop_bhf_TDlearn(ptID) fits BHF data to Q-learning model.
%
%	The output data structure TDdata includes timecourses of Q-learning parameters
%	over trials and fit parameters.
%
%   'whichTrialsLM' input argument specifies, as a string, which trials to
%       evaluate the linear model correlating high gamma activity and TD
%       learning model variables. (default: '40:end')
% 
%   'whichTrialsTD' input argument specifies, as a string, which trials to
%       include in the TD learning model. 
%       options are: 'red', 'orange', 'yellow', and 'control' (default: 'all')
%
%   The fourth input argument is a boolean, and specifies whether to
%   plot all significant electrodes. (default: false)
% score
set(0,'defaultfigurerender','painters'); % making sure we can edit figures after saving.

if nargin<3
    whichTrialsLM = '1:end';
    whichTrialsTD = 'all';
    plotFlag = false;
end

% author: EHS20181005

% ptID = '202215';
% whichRL = 'vanilla';
% whichTrialsLM = '1:end';
% plotFlag = true;

parentDir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev'];
nevList = dir(parentDir)
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,~,~,anatomicalLocs,~] = ptTrodesBART_2(ptID);

% loading behavioral matFile
matFile = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\' ptID '.bartBHV.mat'];
load(matFile)

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,~] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
nChans = length(trodeLabels(isECoG));
Fs = NSX.MetaTags.SamplingFreq;

% timing parameters.
pre = 2;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% task parameters in chronological order..
% respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
% isCTRL = logical([data.is_control]);
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
% inflateTimes = trigTimes(trigs==23);
% pointsPerTrial = [data.points];
if length(balloonTimes)>length(outcomeTimes)
    balloonTimes(end) = [];
end
% if length(inflateTimes)>length(outcomeTimes)
%     inflateTimes(end) = [];
% end
% if length(pointsPerTrial)>length(outcomeTimes)
%     pointsPerTrial(end) = [];
% end

% picking trials to TD model 
if strcmp(whichTrialsTD,'all')
    loopTrials = 1:length(outcomeType);
elseif strcmp(whichTrialsTD,'red')
    loopTrials = find(balloonIDs==1);
elseif strcmp(whichTrialsTD,'orange')
    loopTrials = find(balloonIDs==2);
elseif strcmp(whichTrialsTD,'yellow')
    loopTrials = find(balloonIDs==3);
elseif strcmp(whichTrialsTD,'control')
    loopTrials = find(balloonIDs>10);
end

% trimming the trial numbers if task stopped after a cue. 
if loopTrials(end)>length(outcomeTimes)
    nTrials = length(loopTrials)-1;
else
    nTrials = length(loopTrials);
end

% HG filter
[b,a] = butter(4,[70 160]/(Fs/2));

% epoching data - makes a smaller array than number of trials for
% whichtrialsTD~='all'
for tt = nTrials:-1:1
    % epoch the data here [channels X samples X trials]
    % outcome-aligned
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(loopTrials(tt)))-Fs*pre:floor(Fs*outcomeTimes(loopTrials(tt)))+Fs*post);
    % balloon-aligned
    LFPmat1(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(loopTrials(tt)))-Fs*pre:floor(Fs*balloonTimes(loopTrials(tt)))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        % outcome-aligned
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
        % balloon-aligned
        HGmat1(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat1(ch,:,tt)))));
    end
end

% how much smoothing in time? (50 - 100  ms is usually good for BHF)
% Fs = 1000 samples per second
smoothType = 'movmean';
smoothFactor = Fs./5;

% smoothing high gammma tensors.
HGmat = smoothdata(HGmat,2,smoothType,smoothFactor);
HGmat1 = smoothdata(HGmat1,2,smoothType,smoothFactor);

% saving vars
rewardAndRiskHGs.outcomeHG = HGmat;
rewardAndRiskHGs.cueHG = HGmat1;
rewardAndRiskHGs.tSec = tSec;
rewardAndRiskHGs.timeWin = [-pre post];

% % Summary statistics for each channel in time.
% nBanks = sum(outcomeType==1);
% nPops = sum(outcomeType==2);

% defining model data
baselineNorm = true;
if baselineNorm
    bP = [-1.2 -0.2]; % one second, starting a second and a half before the outcome
    rewardData = squeeze(mean(HGmat(:,tSec > 0.25 & tSec < 1.25,:),2))./squeeze(mean(HGmat(:,tSec > bP(1) & tSec < bP(2),:),2));
    riskData = squeeze(mean(HGmat1(:,tSec > 0.25 & tSec < 1.25,:),2))./squeeze(mean(HGmat1(:,tSec > bP(1) & tSec < bP(2),:),2));
else
    rewardData = squeeze(mean(HGmat(:,tSec > 0.25 & tSec < 1.25,:),2));
    riskData = squeeze(mean(HGmat1(:,tSec > 0.25 & tSec < 1.25,:),2));
end

% fitting the neural data to models.
TDdata = TDlearn(ptID,whichRL,rewardData,riskData,2,whichTrialsLM,loopTrials); % ,fileparts(parentDir)

% finding the best alphas with MLE
[~,bestAlphaReward] = max(TDdata.inverseTemperature);
[~,bestAlphaRisk] = max(TDdata.inverseTemperatureRisk);

% loopng over channels to get indices for significant contacts.
crit = 0.05;
for cc = length(TDdata.neuralFit):-1:1
    % significance.
    sigReward(cc) = TDdata.neuralFit(cc).rewardModel.Coefficients{2,6}<crit;
    sigRewardPE(cc) = TDdata.neuralFit(cc).rewardModel.Coefficients{3,6}<crit;
    sigRisk(cc) = TDdata.neuralFit(cc).riskModel.Coefficients{2,6}<crit;
    sigRiskPE(cc) = TDdata.neuralFit(cc).riskModel.Coefficients{3,6}<crit;

    % just using this loop to put the anatomical labels into TDdata
    try
        TDdata.neuralFit(cc).NMManatomicalLabel = deblank(anatomicalLocs{cc});
    catch
        TDdata.neuralFit(cc).NMManatomicalLabel = 'could not get labels for this electrode';
    end
end

%  decide how many quantiles to look at
% [20220607EHS] (just starting with quartiles)
ps = [0.25 0.5 0.75];
ValueQs = quantile(TDdata.V(bestAlphaReward,:),ps);
RewardPEQs = quantile(TDdata.RewardPE(bestAlphaReward,:),ps);
RiskQs = quantile(TDdata.Vrisk(bestAlphaRisk,:),ps);
RiskPEQs = quantile(TDdata.RiskPE(bestAlphaRisk,:),ps);

if plotFlag
    % visual params
   qCols = viridis(length(ps)+1)./max(max(viridis(length(ps)+1)));
    plotALPHA = 0.3;
    dotSize = 10;

    %% [20220607] plotting data from contacts that are significantly correlated with TD learning model variables
    
    Unsigned = true; 

    if baselineNorm && ~Unsigned
        saveDir  = 'D:\Data\Rhiannon\BART_RLDM_outputs\TDlearn\significantContacts_baselined';
    elseif baselineNorm && Unsigned
       saveDir  = 'D:\Data\Rhiannon\BART_RLDM_outputs\TDlearn\Unsigned_significantContacts_baselined';
    else
        saveDir  = 'D:\Data\Rhiannon\BART_RLDM_outputs\TDlearn\significantContacts';
    end

    for ch2 = 1:nChans
        % plotting responses for all variables, iff theres sig
        if (sigReward(ch2) || sigRewardPE(ch2) || sigRisk(ch2) || sigRiskPE(ch2))

            % setting up figure
            if ishandle(ch2); close(ch2); end
            figure(ch2)

            % plot significant Reward estimate example
            subplot(4,2,3)
            hold on
            for q1 = 1:length(ps)+1
                if q1==1
                    nT = sum(TDdata.V(bestAlphaReward,:)<ValueQs(q1));
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)<ValueQs(q1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)<ValueQs(q1)),[],3));
                elseif q1==(length(ps)+1)
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),[],3));
                    nT = sum(TDdata.V(bestAlphaReward,:)>ValueQs(q1-1));
                else
                    nT = sum(TDdata.V(bestAlphaReward,:)<ValueQs(q1) & TDdata.V(bestAlphaReward,:)>ValueQs(q1-1));
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)<ValueQs(q1) & TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.V(bestAlphaReward,:)<ValueQs(q1) & TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),[],3));
                end
                patch([tSec fliplr(tSec)],[qDataBar+(qDataErr./sqrt(nT)) fliplr(qDataBar-(qDataErr./sqrt(nT)))],qCols(q1,:),'facealpha',plotALPHA,'edgecolor','none')
                plot(tSec,qDataBar,'color',qCols(q1,:))
            end
            hold off

            % plot one deets
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to outcome (s)')
            ylabel('BHF power quantiles')
            title(sprintf('value:t(%d) = %.2f; p = %.2f',TDdata.neuralFit(ch2).rewardModel.Coefficients{2,5},TDdata.neuralFit(ch2).rewardModel.Coefficients{2,4},TDdata.neuralFit(ch2).rewardModel.Coefficients{2,6}));


            % plot signficant reward PE example
            subplot(4,2,7)
            hold on
            for q1 = 1:length(ps)+1
                if q1==1
                    nT = sum(TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1));
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1)),[],3));
                elseif q1==(length(ps)+1)
                    nT = sum(TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1));
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),[],3));
                else
                    nT = sum(TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1) & TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1));
                    qDataBar = squeeze(mean(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1) & TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat(ch2,:,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1) & TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),[],3));
                end
                patch([tSec fliplr(tSec)],[qDataBar+(qDataErr./sqrt(nT)) fliplr(qDataBar-(qDataErr./sqrt(nT)))],qCols(q1,:),'facealpha',plotALPHA,'edgecolor','none')
                plot(tSec,qDataBar,'color',qCols(q1,:))
            end
            hold off

            % plot one deets
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to outcome (s)')
            ylabel('BHF power quantiles')
            title(sprintf('reward PE:t(%d) = %.2f; p = %.2f',TDdata.neuralFit(ch2).rewardModel.Coefficients{3,5},TDdata.neuralFit(ch2).rewardModel.Coefficients{3,4},TDdata.neuralFit(ch2).rewardModel.Coefficients{3,6}));


            % plot significant risk estimate example.
            subplot(4,2,4)
            hold on
            for q1 = 1:length(ps)+1
                if q1==1
                    nT = sum(TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1)),[],3));
                elseif q1==(length(ps)+1)
                    nT = sum(TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),[],3));
                else
                    nT = sum(TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1) & TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1) & TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1) & TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),[],3));
                end
                patch([tSec fliplr(tSec)],[qDataBar+(qDataErr./sqrt(nT)) fliplr(qDataBar-(qDataErr./sqrt(nT)))],qCols(q1,:),'facealpha',plotALPHA,'edgecolor','none')
                plot(tSec,qDataBar,'color',qCols(q1,:))
            end
            hold off

            % plot one deets
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to balloon (s)')
            ylabel('BHF power quantiles')
            title(sprintf('risk: t(%d) = %.2f; p = %.2f',TDdata.neuralFit(ch2).riskModel.Coefficients{2,5},TDdata.neuralFit(ch2).riskModel.Coefficients{2,4},TDdata.neuralFit(ch2).riskModel.Coefficients{2,6}));


            % plot significant Risk PE example
            subplot(4,2,8)
            hold on
            for q1 = 1:length(ps)+1
                if q1==1
                    nT = sum(TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1)),[],3));
                elseif q1==(length(ps)+1)
                    nT = sum(TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),[],3));
                else
                    nT = sum(TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1) & TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1));
                    qDataBar = squeeze(mean(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1) & TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),3));
                    qDataErr = squeeze(std(HGmat1(ch2,:,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1) & TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),[],3));
                end
                patch([tSec fliplr(tSec)],[qDataBar+(qDataErr./sqrt(nT)) fliplr(qDataBar-(qDataErr./sqrt(nT)))],qCols(q1,:),'facealpha',plotALPHA,'edgecolor','none')
                plot(tSec,qDataBar,'color',qCols(q1,:))
            end
            hold off

            % plot one deets
            axis tight square
            xlim([-pre+1 post-1])
            xlabel('time relative to balloon (s)')
            ylabel('BHF power quantiles')
            title(sprintf('risk PE: t(%d) = %.2f; p = %.2f',TDdata.neuralFit(ch2).riskModel.Coefficients{3,5},TDdata.neuralFit(ch2).riskModel.Coefficients{3,4},TDdata.neuralFit(ch2).riskModel.Coefficients{3,6}));


            % plot the VALUE quantiles on the behavioral data to illustrate which trials are averaged over.
            trials = 1:nTrials;
            subplot(4,2,1)
            hold on
            plot(TDdata.V(bestAlphaReward,:),'k')
            for q1 = 1:length(ps)+1
                if q1==1
                    scatter(trials(TDdata.V(bestAlphaReward,:)<ValueQs(q1)),TDdata.V(bestAlphaReward,TDdata.V(bestAlphaReward,:)<ValueQs(q1)),dotSize,qCols(q1,:),'filled')
                elseif q1==(length(ps)+1)
                    scatter(trials(TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),TDdata.V(bestAlphaReward,TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),dotSize,qCols(q1,:),'filled')
                else
                    scatter(trials(TDdata.V(bestAlphaReward,:)<ValueQs(q1) & TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),...
                        TDdata.V(bestAlphaReward,TDdata.V(bestAlphaReward,:)<ValueQs(q1) & TDdata.V(bestAlphaReward,:)>ValueQs(q1-1)),dotSize,qCols(q1,:),'filled')
                end
            end
            hold off
            % axis details.
            axis tight
            xlabel('trials')
            ylabel('reward value estimate')

            subplot(4,2,2)
            hold on
            plot(TDdata.Vrisk(bestAlphaRisk,:),'k')
            for q1 = 1:length(ps)+1
                if q1==1
                    scatter(trials(TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1)),TDdata.Vrisk(bestAlphaRisk,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1)),dotSize,qCols(q1,:),'filled')
                elseif q1==(length(ps)+1)
                    scatter(trials(TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),TDdata.Vrisk(bestAlphaRisk,TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),dotSize,qCols(q1,:),'filled')
                else
                    scatter(trials(TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1) & TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),...
                        TDdata.Vrisk(bestAlphaRisk,TDdata.Vrisk(bestAlphaRisk,:)<RiskQs(q1) & TDdata.Vrisk(bestAlphaRisk,:)>RiskQs(q1-1)),dotSize,qCols(q1,:),'filled')
                end
            end
            hold off
            % axis details.
            axis tight 
            xlabel('trials')
            ylabel('risk estimate')

            % title of first plot is
            title([deblank(trodeLabels{ch2}) ' -- ' TDdata.neuralFit(ch2).NMManatomicalLabel])

            % plot the RPE quantiles on the behavioral data to illustrate which trials are averaged over.
            trials = 1:nTrials;
            subplot(4,2,5)
            hold on
            plot(TDdata.RewardPE(bestAlphaReward,:),'k')
            for q1 = 1:length(ps)+1
                if q1==1
                    scatter(trials(TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1)),TDdata.RewardPE(bestAlphaReward,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1)),dotSize,qCols(q1,:),'filled')
                elseif q1==(length(ps)+1)
                    scatter(trials(TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),TDdata.RewardPE(bestAlphaReward,TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),dotSize,qCols(q1,:),'filled')
                else
                    scatter(trials(TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1) & TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),...
                        TDdata.RewardPE(bestAlphaReward,TDdata.RewardPE(bestAlphaReward,:)<RewardPEQs(q1) & TDdata.RewardPE(bestAlphaReward,:)>RewardPEQs(q1-1)),dotSize,qCols(q1,:),'filled')
                end
            end
            hold off
            % axis details.
            axis tight
            xlabel('trials')
            ylabel('reward prediction error estimate')

            subplot(4,2,6)
            hold on
            plot(TDdata.RiskPE(bestAlphaRisk,:),'k')
            for q1 = 1:length(ps)+1
                if q1==1
                    scatter(trials(TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1)),TDdata.RiskPE(bestAlphaRisk,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1)),dotSize,qCols(q1,:),'filled')
                elseif q1==(length(ps)+1)
                    scatter(trials(TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),TDdata.RiskPE(bestAlphaRisk,TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),dotSize,qCols(q1,:),'filled')
                else
                    scatter(trials(TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1) & TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),...
                        TDdata.RiskPE(bestAlphaRisk,TDdata.RiskPE(bestAlphaRisk,:)<RiskPEQs(q1) & TDdata.RiskPE(bestAlphaRisk,:)>RiskPEQs(q1-1)),dotSize,qCols(q1,:),'filled')
                end
            end
            hold off
            % axis details.
            axis tight
            xlabel('trials')
            ylabel('risk prediction error estimate')

            % saving and closing.
            halfMaximize(ch2,'page')
            saveas(ch2,fullfile(saveDir,[ptID '_' deblank(trodeLabels{ch2}) '_exampleResponses_' smoothType '.pdf']))
            close(ch2)
        end
    end % looping over channels
end % if plotting

end % eof


