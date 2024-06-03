function [BARTstats] = BART_colorTracking_spectra(ptID)
% BART_COLORTRACKING_SPECTRA analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_colorTracking_spectra(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005

% patient
if ~exist('ptID','var')
    ptID = '201902r';
end

% get nev
nevList = dir(sprintf('~/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
selectedChans = find(isECoG);
nChans = length(selectedChans);
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% de-noising and resampling LFP at Fnew sampling frequency
notchFilter = true;
removePC = false;
Fnew = 500;
for ch = 1:nChans
    if notchFilter
        [b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        tmp = filtfilt(b,a,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
        data2K(ch,:) = filtfilt(b2,a2,tmp);
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end
clear NSX NEV
Fs = Fnew;
if removePC
    data2K = remove1stPC(data2K);
end

% timing parameters.
pre = 4;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpass = [1 250]; % pick a value that's 50 higher for wavelets
params.Fs = Fs;
params.dBconversion = false;
params.normalized = false;
params.theoreticalNorm = false;
params.baseline = true;

% picking pre-bird baseline period
if params.baseline
    baselineType = 'preCue';
    switch baselineType
        case {'preCue'}
            bP = [-3.5 -2.5];
        case {'preTask'}
            % normalizing based on the mean spectrum from a pre-task baseline epoch.
            secsPreTask = 50;
            if (ppData.Event.trigTimes(1)./3e4)<secsPreTask;
                display(sprintf('only %d seconds before first trigger.',ppData.Event.trigTimes(1)./3e4));
                % TODO:: then do spectrum
            else
                % TODO:: do spectrum.
            end
    end
end

% task parameters in chronological order..
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% task parameters in chronological order..
respTimes = trigTimes(trigs==23);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);

% adjusting balloon index to match alignment index
balloonType = balloonType(1:nTrials);

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting bank/pop responses
for ch2 = 1:nChans
    % calculate baseline spectra
    
    % epoch the spectral data for each channel.\
    for tt = 1:nTrials
        updateUser('finished spectral calculations',tt,50,nTrials);
        % epoch the data here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
        
        % do spectral claculations here
        [W,period,scale] = basewaveERP(LFPmat(ch2,:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
        
        if params.dBconversion
            normalization = 'decibels';
            Sft(:,:,tt) = abs(10*log10((W)));
        elseif params.baseline
            Sft(:,:,tt) = abs(W);
        elseif params.normalized
            normalization = 'normByFreq';
            Sft(:,:,tt) = normlogspec(abs((W))')';
        elseif params.theoreticalNorm
            normalization = 'theoreticalNorm';
            Sft(:,:,tt) = abs((W))...
                ./repmat(period',1,size(W,2));
        end
        
        % phase angle.
        Pft(:,:,tt) = angle(W);
        
    end
    
    % doing baseline normalization
    if (params.baseline & baselineType=='preCue')
        normalization = [baselineType 'baselineNorm'];
        tmp = Sft./repmat(mean(mean(Sft(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft,2),size(Sft,3));
        Sft = tmp;
        clear tmp
    elseif (params.baseline & baselineType=='preTask')
        normalization = [baselineType 'baselineNorm'];
        % TODO:: add pre task baseline normalization.
    end
    
    % determining spectral scales in Hz
    scaleFreqs = period.^-1;
    
    
    % Summary statistics for each channel
    nBanks = sum(outcomeType==1);
    nPops = sum(outcomeType==2);
    
    %     % permutation tests between conditions for each channel [NOT CLUSTER
    %     % CORRECTING OVER CHANNELS>>>]
    %     nPerms = 1000;
    %     for s = 1:size(HGmat,2)
    %         updateUser(sprintf('channel %d: running permutation test for time sample',ch2),s,100,size(HGmat,2))
    %         p(ch2,s) = permtest2(HGmat(ch2,s,outcomeType==1),HGmat(ch2,s,outcomeType==2),0,nPerms);
    %     end
    
    % visualize the data here for each channel
    if ishandle(ch2); close(ch2); end
    fig1=figure(ch2);
    fig1.Renderer='Painters'; % this line is for vectorized graphics
    
    
    % power colormap
    colormap(parula)
    
    % x axis limits.
    xlims = [-pre+1 post-1];
    
    % ERPs ## NOTE THAT CONTROLS AREN"T INCLUDED AT THIS POINT ##
    subplot(3,2,1)
    hold on
    plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(LFPmat(ch2,tSec>=-pre+1 & tSec<=post-1,balloonType==1),3),'color',cMap(1,:))
    plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(LFPmat(ch2,tSec>=-pre+1 & tSec<=post-1,balloonType==2),3),'color',cMap(2,:))
    plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(LFPmat(ch2,tSec>=-pre+1 & tSec<=post-1,balloonType==3),3),'color',cMap(3,:))
    plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(LFPmat(ch2,tSec>=-pre+1 & tSec<=post-1,balloonType==4),3),'color',cMap(4,:))
    
    % labels
    axis square
    xl = xlim;
    yl = ylim;
    %     text(xl(2)-1,yl(2)-10,'red','color',rgb('forestgreen'),'fontweight','bold')
    %     text(xl(2)-1,yl(2)-15,'orange','color',rgb('orangered'),'fontweight','bold')
    %     text(xl(2)-1,yl(2)-20,'yellow','color',rgb('forestgreen'),'fontweight','bold')
    %     text(xl(2)-1,yl(2)-25,'gray','color',rgb('forestgreen'),'fontweight','bold')
    hold off
    xlim(xlims)
    xlabel('time relative to outcome (s)')
    ylabel('ERP (uV)')
    
    % beta from the spectrogram.
    subplot(3,2,3)
    hold on
    for cnd = 1:4
        patch([tSec(tSec>=-pre+1 & tSec<=post-1),fliplr(tSec(tSec>=-pre+1 & tSec<=post-1))],...
            cat(2,mean(mean(Sft(scaleFreqs>15 & scaleFreqs<25,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd)),3)...
            +(std(mean(Sft(scaleFreqs>15 & scaleFreqs<25,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd)),[],3)./sum(balloonType==cnd | balloonType==10+cnd)),...
            fliplr(mean(mean(Sft(scaleFreqs>15 & scaleFreqs<25,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd)),3)...
            -(std(mean(Sft(scaleFreqs>15 & scaleFreqs<25,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd)),[],3)./sum(balloonType==cnd | balloonType==10+cnd)))),...
            cMap(cnd,:),'edgecolor','none','facealpha',0.3)
        plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(mean(Sft(scaleFreqs>15 & scaleFreqs<25,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),3),'color',cMap(cnd,:),'linewidth',1)
    end
    % labels
    axis square tight
    xl = xlim;
    yl = ylim;
    %     text(xl(2)-1,yl(2)-0.2,'banked','color',rgb('forestgreen'),'fontweight','bold')
    %     text(xl(2)-1,yl(2)-0.4,'popped','color',rgb('orangered'),'fontweight','bold')
    hold off
    xlim(xlims)
    xlabel('time relative to outcome (s)')
    ylabel('beta power (re: baseline)')
    
    % teta from the spectrogram.
    subplot(3,2,5)
    hold on
    for cnd = 1:4
        patch([tSec(tSec>=-pre+1 & tSec<=post-1) fliplr(tSec(tSec>=-pre+1 & tSec<=post-1))],...
            cat(2,mean(mean(Sft(scaleFreqs>4 & scaleFreqs<8,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),3)...
            +(std(mean(Sft(scaleFreqs>4 & scaleFreqs<8,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),[],3)./sum(balloonType==cnd | balloonType==10+cnd)),...
            fliplr(mean(mean(Sft(scaleFreqs>4 & scaleFreqs<8,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),3)...
            -(std(mean(Sft(scaleFreqs>4 & scaleFreqs<8,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),[],3)./sum(balloonType==cnd | balloonType==10+cnd)))),...
            cMap(cnd,:),'edgecolor','none','facealpha',0.3)
        plot(tSec(tSec>=-pre+1 & tSec<=post-1),mean(mean(Sft(scaleFreqs>4 & scaleFreqs<8,tSec>=-pre+1 & tSec<=post-1,balloonType==cnd | balloonType==10+cnd),1),3),'color',cMap(cnd,:),'linewidth',1)
    end
    % labels
    axis square tight
    xl = xlim;
    yl = ylim;
    %     text(xl(2)-1,yl(2)-0.2,'banked','color',rgb('forestgreen'),'fontweight','bold')
    %     text(xl(2)-1,yl(2)-0.4,'popped','color',rgb('orangered'),'fontweight','bold')
    hold off
    xlim(xlims)
    xlabel('time relative to outcome (s)')
    ylabel('theta  power (re: baseline)')
    
    % spectra
    for cnd = 1:3
        ax3 = subplot(3,2,2*cnd);
        surf(tSec(tSec>=-pre+1 & tSec<=post-1),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=-pre+1 & tSec<=post-1),balloonType==cnd | balloonType==10+cnd),3))),'edgecolor','none');
        set(gca,'Yscale','log');
        view(2)
        axis xy tight square
        colormap(ax3,parula)
        h = colorbar;
        ylabel(h, 'power re: baseline')
        xlim(xlims)
        ylim([0 params.fpass(2)-50])
        xlabel('time relative to outcome (s)')
        ylabel('frequency (Hz)')
        switch cnd
            case 1
                title([deblank(trodeLabels{ch2}) ' [yellow balloons]'])
            case 2
                title([deblank(trodeLabels{ch2}) ' [orange balloons]'])
            case 3
                title([deblank(trodeLabels{ch2}) ' [red balloons]'])
        end
    end
    
    
    %     % gray spectra
    %     ax3 = subplot(3,2,5);
    %     surf(tSec(tSec>=-pre+1 & tSec<=post-1),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=-pre+1 & tSec<=post-1),balloonType==14),3))),'edgecolor','none');
    %     set(gca,'Yscale','log');
    %     view(2)
    %     axis xy tight square
    %     colormap(ax3,parula)
    %     h = colorbar;
    %     ylabel(h, 'power re: baseline')
    %     xlim(xlims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to outcome (s)')
    %     ylabel('frequency (Hz)')
    %     title([deblank(trodeLabels{ch2}) ' [gray balloons]'])
    %
    
    
    
    %     % plotting phases.
    %     % loading phase colormap
    %     load('/home/user1/code/matlab/analyzeMSIT/dependencies/phaseColormap.mat')
    %
    %     % bank phases
    %     ax5 = subplot(3,2,5);
    %     surf(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==1),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==1),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax5,phaseCmap)
    %     xlim(xlims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to outcome (s)')
    %     ylabel('frequency (Hz)')
    %     title([deblank(trodeLabels{ch2}) ' [points banked]'])
    %
    %     % pop phases
    %     ax6 = subplot(3,2,6);
    %     surf(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==2),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==2),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax6,phaseCmap)
    %     xlim(xlims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to outcome (s)')
    %     title([deblank(trodeLabels{ch2}) ' [balloon pops]'])
    
    % overall deets
    halfMaximize(ch2,'left')
    
    % saving data and figures
    saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s/',ptID);
    if ~exist(saveDir,'dir')
        [~,msg] = mkdir(saveDir)
    end
    % print([saveDir '/' sprintf('pt%s_%s_banksvspops_spectrograms_%s',ptID,deblank(trodeLabels{ch2}),normalization)],'-dpdf','-fillpage')
    saveas(ch2,[saveDir '/' sprintf('pt%s_%s_colorTracking_spectrograms_%s',ptID,deblank(trodeLabels{ch2}),normalization) '.pdf']);
    
    close(ch2)
    
end


% do statistics for each channel

% visualize the data here for each channel







