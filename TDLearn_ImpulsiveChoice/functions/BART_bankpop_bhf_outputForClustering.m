function [rewardData,riskData] = BART_bankpop_bhf_outputForClustering(ptID)   % BARTstats
% BART_BANKPOP_BHF_OUTPUTFORCLUSTERING outputs cue and outcome aligned high
%   gamma data for clustering. 
%

% author: EHS20220922
% ptID = '202201';


parentDir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev'];
nevList = dir(parentDir)
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% loading behavioral matFile
matFile = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\' ptID '.bartBHV.mat'];
load(matFile)

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
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

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = logical([data.is_control]);
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
pointsPerTrial = [data.points];
if length(balloonTimes)>length(outcomeTimes)
    balloonTimes(end) = [];
end
if length(inflateTimes)>length(outcomeTimes)
    inflateTimes(end) = [];
end
if length(pointsPerTrial)>length(outcomeTimes)
    pointsPerTrial(end) = [];
end
nTrials = length(outcomeType);

% HG filter
[b,a] = butter(4,[70 160]/(Fs/2));

% epoching data
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    % outcome-aligned
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    % balloon-aligned
    LFPmat1(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
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

% Summary statistics for each channel in time.
nBanks = sum(outcomeType==1);
nPops = sum(outcomeType==2);

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
TDdata = TDlearn(ptID,whichRL,rewardData,riskData,2,fileparts(parentDir));

