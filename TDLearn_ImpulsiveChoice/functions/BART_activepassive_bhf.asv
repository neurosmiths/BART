
function [] = NBACK_bhf(ptID)   % NBACKstats
% NBACK_nback_bhf analyzes and visualizes LFP for the NBACK task.
%
%   [NBACKstats] = NBACK_bhf(ptID) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20191028
% 
% ptID = '201902r';

nevList = dir(['D:\Data\preProcessed\NBACK_preprocessed\' ptID '\Data\' ptID '_N_Back_*.nev']);
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesNBACK(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns6']));

% data parameters
nChans = length(trodeLabels(isECoG));
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop

% number of trials & adjustment if no final outcome occurs. 
nTrials = length(outcomeTimes);
if lt(nTrials,length(balloonTimes))
    balloonType(end) = [];

end
if lt(nTrials,length(inflateTimes))
    inflateTimes(end) = [];
end

% timing parameters.
pre = 3;
post = ceil(median(outcomeTimes-inflateTimes));
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% HG filter
[b,a] = butter(4,[70 200]/(Fs/2));

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
end

% how much smoothing in time? (50 - 100  ms is usually good for BHF) 
smoothFactor = Fs./10;
% Summary statistics for each channel in time. 
nActives = sum(balloonType<4); 
nPassives = sum(balloonType>10);
meanActive = smoothdata(squeeze(mean(HGmat(:,:,balloonType<4),3)),2,'movmean',smoothFactor)';
errActive = smoothdata(squeeze(std(HGmat(:,:,balloonType<4),0,3)),2,'movmean',smoothFactor)'./sqrt(nActives);
meanPassive = smoothdata(squeeze(mean(HGmat(:,:,balloonType>10),3)),2,'movmean',smoothFactor)';
errPassive = smoothdata(squeeze(std(HGmat(:,:,balloonType>10),0,3)),2,'movmean',smoothFactor)'./sqrt(nPassives);


% % permutation tests between conditions for each channel 
% if ~exist(sprintf('%s/BHFclusterStats_activePassive_%s.mat',nevList.folder,ptID),'file')
%     tic
%     [pval,t_orig,clust_info,~,est_alpha] = clust_perm2(HGmat(:,:,balloonType<4),HGmat(:,:,balloonType>10),adjacentChanMat);
%     A = toc; 
%     fprintf('\ncluster-based permutation stats took %d seconds',A)
%     save(sprintf('%s/BHFclusterStats_activePassive_%s.mat',nevList.folder,ptID),'pval','t_orig','clust_info','est_alpha')
% else
%     fprintf('stats have already been done. Aren"t you lucky!')
%     load(sprintf('%s/BHFclusterStats_activePassive_%s.mat',nevList.folder,ptID),'pval','t_orig','clust_info','est_alpha')
% end

% for visualizing the data.  % plotting bank/pop responses
for ch2 = 1:nChans
    
    % visualize the data here for each channel
    if ishandle(ch2); close(ch2); end
    figure(ch2)
    hold on
    patch([tSec fliplr(tSec)],[meanPassive(:,ch2)'-errPassive(:,ch2)' fliplr(meanPassive(:,ch2)'+errPassive(:,ch2)')],rgb('gray'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,meanPassive(:,ch2)','color',rgb('gray'))
    text(0,mean(meanPassive(:,ch2)),sprintf('N = %d',uint16(nPassives)),'color',rgb('gray'))
    patch([tSec fliplr(tSec)],[meanActive(:,ch2)'-errActive(:,ch2)' fliplr(meanActive(:,ch2)'+errActive(:,ch2)')],rgb('black'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,meanActive(:,ch2)','color',rgb('black'))
    text(0,mean(meanActive(:,ch2)),sprintf('N = %d',uint16(nActives)),'color',rgb('black'))
%     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
    hold off
    axis tight square
    xlim([-(pre-1) post])
    xlabel('time relative to balloon (s)')
    ylabel('broadband high frequency power')
    title(deblank(trodeLabels{ch2}))
    
    
    
    % saving figures
    saveDir = sprintf('~/Figs/BART/%s',ptID);
    if exist(saveDir,'dir')
        %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
        saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_activepassive_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
    end
    
end





