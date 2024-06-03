function [] = BART_bankpop_bhf_trajectories(ptID)   % BARTstats
% BART_BANKPOP_BHF analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_bankpop_bhf_trajectories(ptID,nevFile) examines
%   Neural trajectories for each trial.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.

% author: EHS20181005


ptID = '201910';


nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

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
pre = 3;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);
tMap = colormap(redblue(length(tSec)));

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);

% HG filter
[b,a] = butter(4,[70 200]/(Fs/2));

% how much smoothing in time? (50 - 100  ms is usually good for BHF) ./4
% for gaussian
smoothFactor = Fs./20;
timeIdcs = tSec>-0.5 & tSec<1;

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
    
    % smoothing data with 100 ms kernel.
    tmp = smoothdata(squeeze(HGmat(:,:,tt))','gaussian',smoothFactor)';
    
    % z-scoring HGmat for each trial and doing PCA
    %     [coeff,score,latent,~,explained] = pca(zscore(squeeze(HGmat(:,:,tt)),[],2)');
    
    % not z-scoring
    [coeff(:,:,tt),score(:,:,tt),latent(:,:,tt),~,explained(:,tt)] = pca(tmp');
    
    % plotting each trial trajectory
    plotEachTrial = true;
    if plotEachTrial
        figure(1)
        hold on
        if outcomeType(tt)==1
            plot3(squeeze(score(timeIdcs,1,tt)),squeeze(score(timeIdcs,2,tt)),squeeze(score(timeIdcs,3,tt)),'color',rgb('palegreen'),'linewidth',0.25)
        elseif outcomeType(tt)==2
            plot3(squeeze(score(timeIdcs,1,tt)),squeeze(score(timeIdcs,2,tt)),squeeze(score(timeIdcs,3,tt)),'color',rgb('coral'),'linewidth',0.25)
        end
        hold off
    end
end

% plotting mean trajectories.
plotMeans = true;
if plotMeans
    hold on
    plot3(squeeze(mean(score(timeIdcs,1,outcomeType==1),3)),squeeze(mean(score(timeIdcs,2,outcomeType==1),3)),squeeze(mean(score(timeIdcs,3,outcomeType==1),3)),'color',rgb('forestgreen'),'linewidth',2)
    plot3(squeeze(mean(score(timeIdcs,1,outcomeType==2),3)),squeeze(mean(score(timeIdcs,2,outcomeType==2),3)),squeeze(mean(score(timeIdcs,3,outcomeType==2),3)),'color',rgb('orangered'),'linewidth',2)
    hold off
end

% saving figures
saveDir = sprintf('~/data/BART/BART_EMU/%s/Figs',ptID);
if exist(saveDir,'dir')
    saveas(1,fullfile(saveDir,sprintf('pt%s_banksvspops_BHF_trajectories.pdf',ptID)))
end








