function [TDdata] = TDlearn(ptID,whichTD,rewardData,riskData,startTrial,LMtrials,TDtrials,varargin)
% TDLEARN fits a temporal difference learning RL algorithm to BART choice data.
%
%	This function is useful for fitting neural data to Q-learning curves.
%
%	ptID is a patient identifier that will get BART behavior
%	whichTD is a string that indicates whether to use 'vanilla' or 'distributional' RL models.
%	rewardData is a matrix of neural data or responses for each trial and channel to fit to reward model.
%	riskData is a matrix of neural data or responses for each trial and channel to fit to risk model.
%	startTrial specifies the first ttrial from which iteration begins
%   LMtrials is a matlab formatted string indicating which trials to
%       evaluate a linear model across (e.g. '40:80' or 'end-40:end')
%   TDtrials specifies which trials to evaluate the TD model 
%
%   A sixth optional input argument specifies a folder to save temporal
%   difference learning model fit structures.

% author: EHS20201102
% edited: RC

if nargin<3
    error('function requires data')
elseif nargin<5
    fprintf('\nStarting trial not specified. Looking at correlations starting from trial 2')
    startTrial = 2;
elseif nargin<6
    LMtrials = '40:end';
elseif nargin==8
    saveDir = varargin{1};
end

parentDir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev'];
nevList = dir(parentDir);nevFile = fullfile(nevList.folder,nevList.name);
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART_2(ptID);
% trodeLabels = ptTrodesBART(ptID);

% initializing bhv output
TDdata.patientID = ptID;
TDdata.type = whichTD;

% loading matFile
matFile = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\' ptID '.bartBHV.mat'];
load(matFile)
dataBHV = data;
clear data

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% [20170713] I made a small error in the handle_input2.m script such that
% sometimes there is an additional 23 after the initial 23 (signifying the
% start of the inflation). The next two lines remove that second 23.
infIdx = trigs==23;
infIdx([diff(infIdx)==0; false] & infIdx==1) = 0;

% trial IDs
isCTRL = logical([dataBHV.is_control]);
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24;
outcomeTypeCat(outcomeType==1) = {'A'};
outcomeTypeCat(outcomeType==2) = {'B'};
categorical(outcomeTypeCat);

% response times
respTimes = trigTimes(trigs==25 | trigs==26);
% trialStarts = trigTimes(trigs<15);
% trialStarts = trialStarts(2:length(respTimes)+1);
% inflateStarts = trigTimes(infIdx);
if ~exist('TDtrials','var')
    nTrials = length(respTimes);
    % adjusting for numbers of trials.
    balloonIDs = balloonIDs(1:nTrials);
    isCTRL = isCTRL(1:nTrials);
else
    nTrials = length(TDtrials);
    % adjusting for numbers of trials.
    balloonIDs = balloonIDs(TDtrials);
    isCTRL = isCTRL(TDtrials);
    outcomeType = outcomeType(TDtrials);
    outcomeTypeCat = outcomeTypeCat(TDtrials);
end

% risk colormap
cMap(2,:) = [1 0.9 0];
cMap(3,:) = [1 0.5 0];
cMap(4,:) = [1 0 0];
cMap(1,:) = [0.5 0.5 0.5];

% % risk and reward colormap
% rcMap(1,:) = [0.5 0.5 0.5];
% rcMap(2,:) = [1 0 0];
% rcMap(3,:) = [1 0.5 0];
% rcMap(4,:) = [1 0.9 0];
% rcMap(5,:) = rgb('lightcoral');
% rcMap(6,:) = rgb('rosybrown');
% rcMap(7,:) = rgb('violet');

% colormap per trial.
balloonColorMap = ones(length(balloonIDs),3)*.5;
% riskRewardColorMap = zeros(length(balloonIDs),3);

% populating balloon color map
for x = 1:3
    balloonColorMap(balloonIDs==x,:) = repmat(cMap(x+1,:),sum(balloonIDs==x),1);
end

% making sure the matrix is oriented correctly:: [channels X trials]
[a,b] = size(riskData);
if isequal(b,nTrials)
    nChannels = a;
    tmpTrials = b;
elseif isequal(a,nTrials)
    riskData = riskData';
    rewardData = rewardData';
    nChannels = b;
    tmpTrials = a;
else
    error(['Are you sure the neural data variable is the right size? The script is detecting ' nTrials ' trials, but the data appears to have ' tmpTrials ' trials and ' nChannels ' channels'])
end

% points vector
pointsPerTrial = diff([0 [dataBHV.score]]); %added BHV

% fitting TD learning models.
if strcmp(whichTD,'vanilla')
    % Initalize variables.
    Reward = zeros(1,nTrials);
    Risk = zeros(1,nTrials);
    RewardPE = zeros(1,nTrials);
    XPE = zeros(1,nTrials);
    Reward(1) = dataBHV(1).points;
    Risk(1) = 0.5;								% starting with a coin flip for risk

    expectedReward = zeros(1,nTrials);
    expectedRisk = zeros(1,nTrials);
    a = 0.01:0.01:1;		% alphas;
    %b = 1:1:10;             % potential betas

    for as = 1:length(a)
        for t = startTrial:nTrials
            % reward varriable.
            if strcmp(dataBHV(t).result,'banked')
                Reward(t) = pointsPerTrial(t);			% outcome on current trial.
            else
                Reward(t) = 0;
            end

            % risk variable.
            if isCTRL(t)
                Risk(t) = 0;
            else
                Risk(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
                % current trial risk is defined as P(pop) on previously observed balloons.
            end

            % updating risk and reward PE
            RewardPE(t) = Reward(t) - expectedReward(t-1);
            XPE(t) = Risk(t) - expectedRisk(t-1);

            % updating value
            expectedRisk(t) = expectedRisk(t-1) + a(as)*XPE(t-1);
            expectedReward(t) = expectedReward(t-1) + a(as)*RewardPE(t-1);

            % if t==nTrials; keyboard; end;
        end

        whichLink = 'logit';
        % similar but using probit instead of logit. Logit gives the same
        % answer as above.
        B = glmfit(expectedRisk',categorical(outcomeType),'binomial','link',whichLink);
        Brisk = glmfit(expectedReward',categorical(outcomeType),'binomial','link',whichLink);

        % saving variables in a struct for output
        TDdata.nTrials = nTrials;
        TDdata.a = a;
        % inverse temperature parameters from mnrfit
        TDdata.inverseTemperature(as) = B(2);
        TDdata.inverseTemperatureRisk(as) = Brisk(2);
        % TD vars
        TDdata.R(as,:) = Reward;
        TDdata.V(as,:) = expectedReward;
        TDdata.RewardPE(as,:) = RewardPE;
        TDdata.perTrialRisk(as,:) = Risk;
        TDdata.RiskPE(as,:) = XPE;
        TDdata.Vrisk(as,:) = expectedRisk;

    end

    
    %% finding best alphas    
    % adjusting for inverted inverse temperature values.
    if sign(sum(sign(TDdata.inverseTemperatureRisk(5:45))))<0
        TDdata.inverseTemperatureRisk = -TDdata.inverseTemperatureRisk;
    end
    if sign(sum(sign(TDdata.inverseTemperature(5:45))))<0
        TDdata.inverseTemperature = -TDdata.inverseTemperature;
    end  

    % ..for reward
    [~,maxBidx] = max(TDdata.inverseTemperature);
    bestAlpha = a(maxBidx);

    % ... for risk
    [~,maxBidxR] = max(TDdata.inverseTemperatureRisk);
    bestAlphaRisk = a(maxBidxR);


    %% fitting temporal difference learning models to BEHAVIOR
    % TODO:: put all of these signals in a table and run LMMs as controls.
    Xpe = eval(['TDdata.RewardPE(maxBidx,' LMtrials ')']);
    Ype = eval(['TDdata.RiskPE(maxBidxR,' LMtrials ')']);
    Xv = eval(['TDdata.V(maxBidx,' LMtrials ')']);
    Yv = eval(['TDdata.Vrisk(maxBidxR,' LMtrials ')']);
    
    % fitting linear models between risk and reward.
    PElm = fitlm(Xpe',Ype');	% This data appears to cluster rather than fit a line, so going to plot with balloon color
    Vlm = fitlm(Xv',Yv');

    % correlations between reward and risk estimates and PEs.
    TDdata.rewardVsRiskPElm.PEmodel = PElm;
    TDdata.valueVsRiskEstimateLM.ValueModel = Vlm;


    %% coincatentaing a bunch of tables for LMMs (optimal learning rates only)
    % reward ECoG
    tmpTbl1 = array2table(rewardData','VariableNames',trodeLabels(isECoG)');
    % risk ECOG
    tmpTbl3 = array2table(riskData','VariableNames',trodeLabels(isECoG)');

    % Asymmetric Variables
    positiveRewardPE = TDdata.RewardPE(maxBidx,:)';
    positiveRewardPE(sign(TDdata.RewardPE(maxBidx,:)) ~= 1) = 0;

    negativeRewardPE = TDdata.RewardPE(maxBidx,:)';
    negativeRewardPE(sign(TDdata.RewardPE(maxBidx,:)) ~= -1) = 0;

    positiveRiskPE = TDdata.RiskPE(maxBidx,:)';
    positiveRiskPE(sign(TDdata.RiskPE(maxBidx,:)) ~= 1) = 0;

    negativeRiskPE = TDdata.RiskPE(maxBidx,:)';
    negativeRiskPE(sign(TDdata.RiskPE(maxBidx,:)) ~= -1) = 0;
    
    % TD model variables.
    tmpTbl2 = table(TDdata.V(maxBidx,:)',TDdata.RewardPE(maxBidx,:)',TDdata.Vrisk(maxBidxR,:)',TDdata.RiskPE(maxBidxR,:)',outcomeTypeCat', abs(TDdata.RewardPE(maxBidx,:))', abs(TDdata.RiskPE(maxBidxR,:))',...
        positiveRewardPE,negativeRewardPE, positiveRiskPE, negativeRiskPE,...
        'VariableNames',{'ValueEstimate','RewardPE','RiskEstimate','RiskPE','outcome', 'UnsignedRewardPE', 'UnsignedRiskPE', 'positiveRewardPE', 'negativeRewardPE', 'positiveRiskPE', 'negativeRiskPE'});
    rewardTbl = cat(2,tmpTbl1,tmpTbl2);
    riskTbl = cat(2,tmpTbl3,tmpTbl2);


    %% fitting temporal difference learning models to NEURAL DATA
    % fitting the model for each channel for the optimal alpha.
    for ch = 1:nChannels
        TDdata.neuralFit(ch).bestAlpha = bestAlpha;
        TDdata.neuralFit(ch).bestAlphaRisk = bestAlphaRisk;
        TDdata.neuralFit(ch).trodeLabel = tmpTbl1.Properties.VariableNames{ch};

        try
            TDdata.neuralFit(ch).rewardModel = fitglme(rewardTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ ValueEstimate + RewardPE']);
        catch ME
            TDdata.neuralFit(ch).rewardModel = ME.identifier;
        end

        try
            TDdata.neuralFit(ch).riskModel = fitglme(riskTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ RiskEstimate + RiskPE']);
        catch ME
            TDdata.neuralFit(ch).riskModel = ME.identifier;
        end
        
        % controls for balloon pop. 
        balloonPopControl = true;
        if balloonPopControl
            try
                TDdata.neuralFit(ch).rewardModelPOPCTRL = fitglme(rewardTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ ValueEstimate + RewardPE + (ValueEstimate|outcome) + (RewardPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).rewardModelPOPCTRL = ME.identifier;
            end
    
            try
                TDdata.neuralFit(ch).riskModelPOPCTRL = fitglme(riskTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ RiskEstimate + RiskPE + (RiskEstimate|outcome) + (RiskPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).riskModelPOPCTRL = ME.identifier;
            end
            % unsigned
             try
                TDdata.neuralFit(ch).unsignedRewardModelPOPCTRL = fitglme(rewardTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ ValueEstimate + UnsignedRewardPE + (ValueEstimate|outcome) + (UnsignedRewardPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).unsignedRewardModelPOPCTRL = ME.identifier;
            end
    
            try
                TDdata.neuralFit(ch).unsignedRiskModelPOPCTRL = fitglme(riskTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ RiskEstimate + UnsignedRiskPE + (RiskEstimate|outcome) + (UnsignedRiskPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).unsignedRiskModelPOPCTRL = ME.identifier;
            end
            % asymmetric
             try
                TDdata.neuralFit(ch).asymmetricRewardModelPOPCTRL = fitglme(rewardTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ ValueEstimate + positiveRewardPE + negativeRewardPE + (ValueEstimate|outcome) + (positiveRewardPE|outcome) + (negativeRewardPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).asymmetricRewardModelPOPCTRL = ME.identifier;
            end
    
            try
                TDdata.neuralFit(ch).asymmetricRiskModelPOPCTRL = fitglme(riskTbl,[TDdata.neuralFit(ch).trodeLabel ' ~ RiskEstimate + positiveRiskPE + negativeRiskPE + (RiskEstimate|outcome) + (positiveRiskPE|outcome) + (negativeRiskPE|outcome)']);
            catch ME
                TDdata.neuralFit(ch).asymmetricRiskModelPOPCTRL = ME.identifier;
            end
        end

    end

  
    % saving data.
    if nargin==8
        save(fullfile(saveDir,sprintf('%s_TDlearningModelFits.mat',ptID)),'TDdata','-v7.3') 
    end

elseif strcmp(whichTD,'distributional')
    %% not cler wehther it's even worth doing this. 
%     % expectile or quantile learning rule
%     whichDistRL = 'expectile';
% 
%     % fit model to each channel
%     alphas = 0.01:0.01:1;
%     nAlphas = length(alphas);
% 
%     % Initalize variables.
%     Reward(:,1) = repmat(dataBHV(1).points,nAlphas,1);
%     Risk(:,1) = repmat(0.5,nAlphas,1);
%     expectedReward = zeros(nAlphas,nTrials);
%     expectedRisk = zeros(nAlphas,nTrials);
% 
%     % looping over trials.
%     for t = 2:nTrials
%         if strcmp(data(t).result,'banked')
%             Reward(:,t) = repmat(dataBHV(t).points,nAlphas,1);	% outcome on current trial.
%         else
%             Reward(:,t) = zeros(nAlphas,1);
%         end
% 
%         % risk variable.
%         if isCTRL(t)
%             Risk(:,t) = 0;
%         else
%             Risk(:,t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
%             % current trial risk is defined as P(pop) on previously observed balloons.
%         end
% 
%         % update PEs
%         RewardPE(:,t) = Reward(:,t) - expectedReward(:,t-1);
%         XPE(:,t) = Risk(:,t) - expectedRisk(:,t-1);
% 
%         % updating values for risk.
%         if strcmp(whichDistRL,'quantile')
%             expectedRisk(XPE(:,t)>0,t) = expectedRisk(XPE(:,t)>0,t-1) + alphas(XPE(:,t)>0).*ones(sum(XPE(:,t)>0),1);
%             expectedRisk(XPE(:,t)<=0,t) = expectedRisk(XPE(:,t)<=0,t-1) + alphas(XPE(:,t)<=0).*-ones(sum(XPE(:,t)<=0),1);
%         elseif strcmp(whichDistRL,'expectile')
%             expectedRisk(XPE(:,t)>0,t) = expectedRisk(XPE(:,t)>0,t-1) + alphas(XPE(:,t)>0).*XPE(XPE(:,t)>0,t);
%             expectedRisk(XPE(:,t)<=0,t) = expectedRisk(XPE(:,t)<=0,t-1) + alphas(XPE(:,t)<=0).*XPE(XPE(:,t)<=0,t);
%         end
% 
%         % updating values for reward.
%         if strcmp(whichDistRL,'quantile')
%             expectedReward(RewardPE(:,t)>0,t) = expectedReward(RewardPE(:,t)>0,t-1) + alphas(RewardPE(:,t)>0).*ones(sum(RewardPE(:,t)>0),1);
%             expectedReward(RewardPE(:,t)<=0,t) = expectedReward(RewardPE(:,t)<=0,t-1) + alphas(RewardPE(:,t)<=0).*-ones(sum(RewardPE(:,t)<=0),1);
%         elseif strcmp(whichDistRL,'expectile')
%             expectedReward(RewardPE(:,t)>0,t) = expectedReward(RewardPE(:,t)>0,t-1) + alphas(RewardPE(:,t)>0).*RewardPE(RewardPE(:,t)>0,t);
%             expectedReward(RewardPE(:,t)<=0,t) = expectedReward(RewardPE(:,t)<=0,t-1) + alphas(RewardPE(:,t)<=0).*RewardPE(RewardPE(:,t)<=0,t);
%         end
% 
%     end
% 
%     % plotting
%     for xx = 1
%         plt = false;
%         if plt
% 
%             figure
%             subplot(3,2,1)
%             plot(1:nTrials,RewardPE)
%             xlabel('trials')
%             ylabel('reward PE')
%             title(sprintf('patient: %s',ptID))
% 
%             subplot(3,2,2)
%             plot(1:nTrials,expectedReward)
%             xlabel('trials')
%             ylabel('value estimate')
%             title([whichDistRL ' learning rule'])
% 
%             subplot(3,2,3)
%             plot(1:nTrials,XPE)
%             xlabel('trials')
%             ylabel('Risk PE')
% 
%             subplot(3,2,4)
%             plot(1:nTrials,expectedRisk)
%             xlabel('trials')
%             ylabel('risk estimate')
% 
%             % 		subplot(3,2,5)
%             % 		PElm.plot
%             % 		xlabel('reward PE')
%             % 		ylabel('risk PE')
%             % 		title(sprintf('t(%d) = %.2f, p = %.2f, R^2 = %.2f',PElm.DFE,PElm.Coefficients{2,3},PElm.Coefficients{2,4},PElm.Rsquared.Ordinary))
%             %
%             % 		subplot(3,2,6)
%             % 		Vlm.plot
%             % 		xlabel('Vreward')
%             % 		ylabel('Vrisk')
%             % 		title(sprintf('t(%d) = %.2f, p = %.2f, R^2 = %.2f',Vlm.DFE,Vlm.Coefficients{2,3},Vlm.Coefficients{2,4},Vlm.Rsquared.Ordinary))
% 
%             % saving figures
%             halfMaximize(gcf,'left')
%             saveas(gcf,['D:\Data\Rhiannon\BART_RLDM_outputs\TDlearn\' ptID '_' whichDistRL 'TDlearning.pdf'])
%         end
%     end
% 
%     % saving variables in a struct for output
%     TDdata.nTrials = nTrials;
%     TDdata.a = alphas;
%     TDdata.R = Reward;
%     TDdata.V = expectedReward;
%     TDdata.RewardPE = RewardPE;
%     TDdata.perTrialRisk = Risk;
%     TDdata.RiskPE = XPE;
%     TDdata.Vrisk = expectedRisk;

end





