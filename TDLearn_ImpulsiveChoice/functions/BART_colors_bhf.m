function [BARTstats] = BART_colors_bhf(ptID)
% BART_COLORS_BHF analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_colors_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005

ptID = '201901';

nevList = dir(['\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
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
Fs = NSX.MetaTags.SamplingFreq; %MAGIC NUMBER!!!!!!!

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
outcomeTimes = trigTimes(trigs==25 | trigs==26);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
nTrials = length(balloonTimes);


% timing parameters.
allInflationTimes = outcomeTimes - inflateTimes;
pre = 2;
post = ceil(max(allInflationTimes));
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% how much smoothing in time? (50 - 100  ms is usually good for BHF) ./4
% for gaussian
smoothFactor = Fs./20;
timeIdcs = tSec>-0.5 & tSec<12;

% aligning on balloon appearances to start.
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
    
    % smoothing data with 100 ms kernel.
    tmp = smoothdata(squeeze(HGmat(:,:,tt))','gaussian',smoothFactor)';
    
    % not z-scoring
    [coeff(:,:,tt),score(:,:,tt),latent(:,:,tt),~,explained(:,tt)] = pca(tmp');
    
end

% plotting balloon color trajectories. 
figure(1)
subplot(1,2,1)
hold on
plot3(squeeze(mean(score(timeIdcs,1,balloonType==1 | balloonType==11),3)),squeeze(mean(score(timeIdcs,2,balloonType==1 | balloonType==11),3)),squeeze(mean(score(timeIdcs,3,balloonType==1 | balloonType==11),3)),'color',rgb('red'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==2 | balloonType==12),3)),squeeze(mean(score(timeIdcs,2,balloonType==2 | balloonType==12),3)),squeeze(mean(score(timeIdcs,3,balloonType==2 | balloonType==12),3)),'color',rgb('orangeRed'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==3 | balloonType==13),3)),squeeze(mean(score(timeIdcs,2,balloonType==3 | balloonType==13),3)),squeeze(mean(score(timeIdcs,3,balloonType==3 | balloonType==13),3)),'color',rgb('gold'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==4 | balloonType==14),3)),squeeze(mean(score(timeIdcs,2,balloonType==4 | balloonType==14),3)),squeeze(mean(score(timeIdcs,3,balloonType==4 | balloonType==14),3)),'color',rgb('gray'),'linewidth',2)
hold off


% plotting balloon color trajectories. 
subplot(1,2,2)
hold on
plot3(squeeze(mean(score(timeIdcs,1,balloonType==1),3)),squeeze(mean(score(timeIdcs,2,balloonType==1),3)),squeeze(mean(score(timeIdcs,3,balloonType==1),3)),'color',rgb('red'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==2),3)),squeeze(mean(score(timeIdcs,2,balloonType==2),3)),squeeze(mean(score(timeIdcs,3,balloonType==2),3)),'color',rgb('orangeRed'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==3),3)),squeeze(mean(score(timeIdcs,2,balloonType==3),3)),squeeze(mean(score(timeIdcs,3,balloonType==3),3)),'color',rgb('gold'),'linewidth',2)
plot3(squeeze(mean(score(timeIdcs,1,balloonType==4),3)),squeeze(mean(score(timeIdcs,2,balloonType==4),3)),squeeze(mean(score(timeIdcs,3,balloonType==4),3)),'color',rgb('gray'),'linewidth',2)
hold off


keyboard






% analyses flags
includeControls = true;
compareVsControls = false;

smoothFactor = Fs./20;
% plotting bank/pop responses
for ch2 = 1:nChans
    if includeControls
        % do statistics for each channel
        nReds = sum(balloonType==3 | balloonType==13);
        nOranges = sum(balloonType==2 | balloonType==12);
        nYellows = sum(balloonType==1 | balloonType==11);
        nGrays = sum(balloonType==4 | balloonType==14);
        
        meanReds = smooth(squeeze(mean(HGmat(ch2,:,balloonType==3 | balloonType==13),3)),smoothFactor)';
        errReds = smooth(squeeze(std(HGmat(ch2,:,balloonType==3 | balloonType==13),0,3)),smoothFactor)'./sqrt(nReds);
        meanOranges = smooth(squeeze(mean(HGmat(ch2,:,balloonType==2 | balloonType==12),3)),smoothFactor)';
        errOranges = smooth(squeeze(std(HGmat(ch2,:,balloonType==2 | balloonType==12),0,3)),smoothFactor)'./sqrt(nOranges);
        meanYellows = smooth(squeeze(mean(HGmat(ch2,:,balloonType==1 | balloonType==11),3)),smoothFactor)';
        errYellows = smooth(squeeze(std(HGmat(ch2,:,balloonType==1 | balloonType==11),0,3)),smoothFactor)'./sqrt(nYellows);
        meanGrays = smooth(squeeze(mean(HGmat(ch2,:,balloonType==4 | balloonType==14),3)),smoothFactor)';
        errGrays = smooth(squeeze(std(HGmat(ch2,:,balloonType==4 | balloonType==14),0,3)),smoothFactor)'./sqrt(nGrays);
        
    else
        
    end
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
    xlim([-1 2])
    xlabel('time relative to balloon (s)')
    ylabel('broadband high frequency power')
    title(deblank(trodeLabels{ch2}))
    
    % saving data
    saveDir = (['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data']);
    if exist(saveDir,'dir')
        halfMaximize(ch2,'left')
        if compareVsControls
            saveas(ch2,fullfile(saveDir(['pt%s_%s_activeVsControl_BHF.pdf',ptID,deblank(trodeLabels{ch2})])))
        else
            saveas(ch2,fullfile(saveDir(['pt%s_%s_colors_BHF.pdf',ptID,deblank(trodeLabels{ch2})])))
        end
    end
    close(ch2)
    endz





