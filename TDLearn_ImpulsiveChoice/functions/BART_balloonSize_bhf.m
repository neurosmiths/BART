function [BARTstats] = BART_balloonSize_bhf(ptID,BARTdir)
% BART_COLORS_BHF analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_colors_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   [BARTstats] = BART_colors_bhf(ptID,nevFile,chan2Plot) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.

% author: EHS20181005

ptID = '201901'; % this line for debugging only

% [20190315] analysis flag:: 
% if allTrials is true, BHF for alltrials and all active trials.
% if alltrials is false, BHF for each balloon color will be plotted. 
allTrials = true;

% getting data dir
if nargin<2
    BARTdir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
end
nevList = dir([BARTdir '/*.nev']);

% getting electrode information
[trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);

% loading LMM results for balloon size tracking.
tic
whichModel =  'active_each';
whichAnalysis = 'balloonSizeTrackingGLMresults';
fName = sprintf('%s_%s_%sTrodesModel.mat',ptID,whichAnalysis,whichModel);
load(fullfile(BARTdir,fName))
toc

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
Fs = NSX.MetaTags.SamplingFreq; %MAGIC NUMBER!!!!!!!

% task timings
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% only including complete trials; generally => excluding the last trial.
nTrials = length(respTimes);
balloonType = balloonType(1:nTrials);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
balloonIDs = balloonIDs(1:nTrials);
isCTRL = isCTRL(1:nTrials);

% duration of balloon inflation
inflateDurations = respTimes-inflateTimes;

% timing parameters.
pre = 2;
post = 8;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% aligning on balloon appearances to start.
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*inflateTimes(tt))-Fs*pre:floor(Fs*inflateTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
end


% how much do we smooth the data?
% Fs = 1000 for these experiments.
% so this would smooth by 100 ms moving average.
smoothFactor = Fs./10;


% plotting responses
for ch2 = 1:nChans
    includeControls = true;
    if includeControls
        % numbers for each channel
        nReds = sum(balloonType==3 | balloonType==13);
        nOranges = sum(balloonType==2 | balloonType==12);
        nYellows = sum(balloonType==1 | balloonType==11);
        nGrays = sum(balloonType==4 | balloonType==14);
        
        % means and sems for each channel
        meanReds = smooth(squeeze(mean(HGmat(ch2,:,balloonType==3 | balloonType==13),3)),smoothFactor)';
        errReds = smooth(squeeze(std(HGmat(ch2,:,balloonType==3 | balloonType==13),0,3)),smoothFactor)'./sqrt(nReds);
        meanOranges = smooth(squeeze(mean(HGmat(ch2,:,balloonType==2 | balloonType==12),3)),smoothFactor)';
        errOranges = smooth(squeeze(std(HGmat(ch2,:,balloonType==2 | balloonType==12),0,3)),smoothFactor)'./sqrt(nOranges);
        meanYellows = smooth(squeeze(mean(HGmat(ch2,:,balloonType==1 | balloonType==11),3)),smoothFactor)';
        errYellows = smooth(squeeze(std(HGmat(ch2,:,balloonType==1 | balloonType==11),0,3)),smoothFactor)'./sqrt(nYellows);
        meanGrays = smooth(squeeze(mean(HGmat(ch2,:,balloonType==4 | balloonType==14),3)),smoothFactor)';
        errGrays = smooth(squeeze(std(HGmat(ch2,:,balloonType==4 | balloonType==14),0,3)),smoothFactor)'./sqrt(nGrays);
        
    else
        % numbers for each channel
        nReds = sum(balloonType==3);
        nOranges = sum(balloonType==2);
        nYellows = sum(balloonType==1);
        nGrays = sum(balloonType==4);
        
        % means and sems for each channel
        meanReds = smooth(squeeze(mean(HGmat(ch2,:,balloonType==3),3)),smoothFactor)';
        errReds = smooth(squeeze(std(HGmat(ch2,:,balloonType==3),0,3)),smoothFactor)'./sqrt(nReds);
        meanOranges = smooth(squeeze(mean(HGmat(ch2,:,balloonType==2),3)),smoothFactor)';
        errOranges = smooth(squeeze(std(HGmat(ch2,:,balloonType==2),0,3)),smoothFactor)'./sqrt(nOranges);
        meanYellows = smooth(squeeze(mean(HGmat(ch2,:,balloonType==1),3)),smoothFactor)';
        errYellows = smooth(squeeze(std(HGmat(ch2,:,balloonType==1),0,3)),smoothFactor)'./sqrt(nYellows);
        meanGrays = smooth(squeeze(mean(HGmat(ch2,:,balloonType==4),3)),smoothFactor)';
        errGrays = smooth(squeeze(std(HGmat(ch2,:,balloonType==4),0,3)),smoothFactor)'./sqrt(nGrays);
        
    end
    
    
    if allTrials
        % wCTRLS
        HGbar = smooth(squeeze(mean(HGmat(ch2,:,:),3)),smoothFactor)';
        HGerr = smooth(squeeze(std(HGmat(ch2,:,:),0,3)),smoothFactor)'./sqrt(nTrials);
        % woCTRLS        
        HGbar_noCTRLS = smooth(squeeze(mean(HGmat(ch2,:,~isCTRL),3)),smoothFactor)';
        HGerr_noCTRLS = smooth(squeeze(std(HGmat(ch2,:,~isCTRL),0,3)),smoothFactor)'./sqrt(nTrials-sum(isCTRL));
        
        % plotting
        if ishandle(ch2); close(ch2); end
        figure(ch2)
        hold on
        patch([tSec fliplr(tSec)],[HGbar-HGerr fliplr(HGbar+HGerr)],rgb('black'),'facealpha',0.3,'edgecolor','none')
        plot(tSec,HGbar,'color',rgb('black'))
        text(0,2.5,sprintf('N = %d',uint16(nTrials)),'color',rgb('black'))
        patch([tSec fliplr(tSec)],[HGbar_noCTRLS-HGerr_noCTRLS fliplr(HGbar_noCTRLS+HGerr_noCTRLS)],rgb('brown'),'facealpha',0.3,'edgecolor','none')
        plot(tSec,HGbar_noCTRLS,'color',rgb('brown'))
        text(0,2.5,sprintf('N = %d',uint16(nTrials-sum(isCTRL))),'color',rgb('brown'))
        hold off
        
        % deets
        hold off
        axis tight square
        xlim([-pre+1 post-1])
        xlabel('time relative to balloon (s)')
        ylabel('broadband high frequency power')
        title(sprintf('%s (%s) LMM: p = %.2f, t(%d) = %.2f',...
            deblank(trodeLabels{ch2}),anatomicalLocs{ch2},trodeLMM(ch2).model_p,trodeLMM(ch2).model_df,trodeLMM(ch2).model_t))
        
        % saving data
        saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s',ptID);
        if exist(saveDir,'dir')
            halfMaximize(ch2,'left')
            saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        end
        close(ch2)
        
    else
        % visualize the data here for each channel
        if ishandle(ch2); close(ch2); end
        figure(ch2)
        hold on
        patch([tSec fliplr(tSec)],[meanReds-errReds fliplr(meanReds+errReds)],rgb('red'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanReds,'color',rgb('red'))
        text(0,2.5,sprintf('N = %d',uint16(nReds)),'color',rgb('red'))
        patch([tSec fliplr(tSec)],[meanOranges-errOranges fliplr(meanOranges+errOranges)],rgb('orange'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanOranges,'color',rgb('orange'))
        text(0,2.5,sprintf('N = %d',uint16(nOranges)),'color',rgb('orange'))
        patch([tSec fliplr(tSec)],[meanYellows-errYellows fliplr(meanYellows+errYellows)],rgb('yellow'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanYellows,'color',rgb('yellow'))
        text(0,2.5,sprintf('N = %d',uint16(nReds)),'color',rgb('yellow'))
        patch([tSec fliplr(tSec)],[meanGrays-errGrays fliplr(meanGrays+errGrays)],rgb('gray'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanGrays,'color',rgb('gray'))
        text(0,2.5,sprintf('N = %d',uint16(nGrays)),'color',rgb('gray'))
        
        % deets
        hold off
        axis tight square
        xlim([-pre+1 post-1])
        xlabel('time relative to balloon (s)')
        ylabel('broadband high frequency power')
        if ~isempty(anatomicalLocs)
            title(sprintf('%s (%s) LMM: p = %.2f, t(%d) = %.2f',...
                deblank(trodeLabels{ch2}),anatomicalLocs{ch2},trodeLMM(ch2).model_p,trodeLMM(ch2).model_df,trodeLMM(ch2).model_t))
        else
            title(sprintf('%s LMM: p = %.2f, t(%d) = %.2f',...
                deblank(trodeLabels{ch2}),trodeLMM(ch2).model_p,trodeLMM(ch2).model_df,trodeLMM(ch2).model_t))
        end
        
        % saving data
        saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s',ptID);
        if exist(saveDir,'dir')
            halfMaximize(ch2,'left')
            if includeControls
                saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_colors_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
            else
                saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_colors_BHF_noCTRL.pdf',ptID,deblank(trodeLabels{ch2}))))
            end
        end
        close(ch2)
    end
end





