function [BARTstats] = BART_balloonOnset_spectra(ptID)
% BART_BALLOONONSET_SPECTRA analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_balloonOnset_spectra(ptID,nevFile) analyzes LFP data
%   for the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005

% ptID = '202002';


nevList = dir(['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);

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

% resampling LFP at Fnew sampling frequency
notchFilter = true;
Fnew = 500;
for ch = 1:nChans
    if notchFilter
        [b1,a1] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        tmp(ch,:) = filtfilt(b1,a1,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
        
        [b2,a2] = iirnotch(120/(Fnew/2),(120/(Fnew/2))/50);
        data2K(ch,:) = filtfilt(b2,a2,tmp(ch,:));
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end
clear NSX NEV
Fs = Fnew;

% common average re-referencing.
%   data2K = remove1stPC(data2K);

% timing parameters.
pre = 2;
post = 4;
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
            bP = [-2.5 -1.5];
        case {'preTask'}
            % normalizing based on the mean spectrum from a pre-task baseline epoch.
            secsPreTask = 50;
            if (ppData.Event.trigTimes(1)/3e4)<secsPreTask
                fprintf('only %d seconds before first trigger.',ppData.Event.trigTimes(1)./3e4);
                % TODO:: then do spectrum
            else
                % TODO:: do spectrum.
            end
    end
end

% task parameters in chronological order..
% There aren't any trigs that == 4
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
if length(balloonTimes)>length(inflateTimes)
    balloonTimes(end) = [];
end


% balloon color colormap - [yellow, orange, read, gray]
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% task parameters in chronological order..
respTimes = trigTimes(trigs==23);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);


if length(balloonType)>nTrials; balloonType=balloonType(1:end-1); end


%% TODO:: don't analyz 'NaC' trodes...

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting bank/pop responses
for ch2 = 1:nChans
    % epoch the spectral data for each channel.
    for tt = 1:nTrials
        updateUser('finished spectral calculations',tt,50,nTrials);
        % epoch the data here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    end
end
clear ch2 tt

% % reject trials with epileptiform discharges.
% [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionUtility(ptID,LFPmat,false);
% fprintf('\n%d trials out of %d have interictal discharges (%d percent).'...
%     ,length(trialsAcrossChannels),nTrials,(length(trialsAcrossChannels)./nTrials)*100)



plotRawTrials = false;
if plotRawTrials
    for chz = 45:48 %size(LFPmat,1):-1:1

        IEDtrials = outliers(range(squeeze(LFPmat(chz,:,:))));

        figure(chz)
        hold on
        for trl = 1:size(LFPmat,3)/2
            if IEDtrials(trl)
            plot(tSec,LFPmat(chz,:,trl)+(500*trl),'r')
            else
            plot(tSec,LFPmat(chz,:,trl)+(500*trl),'k')
            end
        end
        hold off
        xlabel('time (s)')
        title('electrode label here')
    end

end





for ch2 = 1:nChans
    % epoch the spectral data for each channel.\
    for tt = 1:nTrials
%         if ~ismember(tt,trialsAcrossChannels)
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
                    ./repmat(1./linspace(params.fpass(1),params.fpass(2),length(scale))',1,size(W,2));
            end
            
            % phase angle.
            Pft(:,:,tt) = angle(W);
%         end
    end
    
%     % now removing the empty bad trials. 
%     Sft(:,:,trialsAcrossChannels) = [];
%     Pft(:,:,trialsAcrossChannels) = [];
    
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
    
    
    %% visualize the data here for each channel
    if ishandle(ch2); close(ch2); end
    fig1=figure(ch2);
    fig1.Renderer='Painters'; % this line if you wnat vectorized graphics
    
    % x axis limits.
    xLims = [-1 3];
    %     zlims =
    
    
    %% Hypothesis testing..
    % 1) most basic: significant difference from baseline
    %   - could do permutation test.
    [erH,erP,erCI,erSTATS] = ttest(squeeze(nanmean(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-1 & tSec<=0,:),1),2)),squeeze(nanmean(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=0.25 & tSec<=1.25,:),1),2)));
    
    % 2) ANOVA to look at significant differences among balloon colors
    X = squeeze(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=0 & tSec<=mean(inflateTimes-balloonTimes),balloonType<4),1));
    group = balloonType(balloonType<4);
    [colorP,Tbl,colorStats] = kruskalwallis(X,group,'off');
    
    % saving statistics
    balloonOnsetStats(ch2).contactlabel = deblank(trodeLabels{ch2});
    balloonOnsetStats(ch2).changeFromBaseline.timeWin = [-1 0; 0.25 1.25];
    balloonOnsetStats(ch2).changeFromBaseline.H = erH;
    balloonOnsetStats(ch2).changeFromBaseline.p = erP;
    balloonOnsetStats(ch2).changeFromBaseline.CI = erCI;
    balloonOnsetStats(ch2).changeFromBaseline.stats = erSTATS;
    balloonOnsetStats(ch2).sigColorDiffs.timeWin = [0 mean(inflateTimes-balloonTimes)];
    balloonOnsetStats(ch2).sigColorDiffs.p = colorP;
    balloonOnsetStats(ch2).sigColorDiffs.H = colorP<0.05;
    balloonOnsetStats(ch2).sigColorDiffs.tbl = Tbl;
    balloonOnsetStats(ch2).sigColorDiffs.stats = colorStats;
    
    
    %% plotting...
    % ERPs
    subplot(3,2,1)
    hold on
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(LFPmat(ch2,tSec>=xLims(1) & tSec<=xLims(2),balloonType==1),3),'color',cMap(1,:))
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(LFPmat(ch2,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),3),'color',cMap(2,:))
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(LFPmat(ch2,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),3),'color',cMap(3,:))
    %     plot(tSec(tSec>=-2 & tSec<=2),nanmean(LFPmat(ch2,tSec>=-2 & tSec<=2,balloonType==14),3),'color',cMap(4,:))
    
    % labels
    axis square
    hold off
    xlim(xLims)
    xlabel('time (s)')
    ylabel('ERP (uV)')
    title(ptID)
    
    % BHF from the spectrogram.
    subplot(3,2,2)
    hold on
    % individual lines
    %     plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-2 & tSec<=2,balloonType==1),1),2,'movmean',2e2)),'color',cMap(1,:),'linewidth',0.4)
    %     plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-2 & tSec<=2,balloonType==2),1),2,'movmean',2e2)),'color',cMap(2,:),'linewidth',0.4)
    %     plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-2 & tSec<=2,balloonType==3),1),2,'movmean',2e2)),'color',cMap(3,:),'linewidth',0.4)
    %     plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-2 & tSec<=2,balloonType==4),1),2,'movmean',2e2)),'color',cMap(4,:),'linewidth',0.4)
    
    % patches, rather than single trials.
    patch([tSec(tSec>=xLims(1) & tSec<=xLims(2)) fliplr(tSec(tSec>=xLims(1) & tSec<=xLims(2)))],...
        [nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==1),1),2,'movmean',1e2),3) +  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==1),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==1)) ...
        fliplr(nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==1),1),2,'movmean',1e2),3) -  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(2) & tSec<=xLims(2),balloonType==1),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==1)))],...
        cMap(1,:),'edgecolor','none','facealpha',0.4)
    patch([tSec(tSec>=xLims(1) & tSec<=xLims(2)) fliplr(tSec(tSec>=xLims(1) & tSec<=xLims(2)))],...
        [nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),1),2,'movmean',1e2),3) +  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==2)) ...
        fliplr(nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),1),2,'movmean',1e2),3) -  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==2)))],...
        cMap(2,:),'edgecolor','none','facealpha',0.4)
    patch([tSec(tSec>=xLims(1) & tSec<=xLims(2)) fliplr(tSec(tSec>=xLims(1) & tSec<=xLims(2)))],...
        [nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),1),2,'movmean',1e2),3) +  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==3)) ...
        fliplr(nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),1),2,'movmean',1e2),3) -  std(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),1),2,'movmean',1e2),[],3)./sqrt(sum(balloonType==3)))],...
        cMap(3,:),'edgecolor','none','facealpha',0.4)
    
    % mean lines
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==1),1),2,'movmean',1e2),3),'color',cMap(1,:),'linewidth',2)
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==2),1),2,'movmean',1e2),3),'color',cMap(2,:),'linewidth',2)
    plot(tSec(tSec>=xLims(1) & tSec<=xLims(2)),nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=xLims(1) & tSec<=xLims(2),balloonType==3),1),2,'movmean',1e2),3),'color',cMap(3,:),'linewidth',2)
    %     plot(tSec(tSec>=-2 & tSec<=2),nanmean(smoothdata(nanmean(Sft(scaleFreqs>70 & scaleFreqs<200,tSec>=-2 & tSec<=2,balloonType==4),1),2,'movmean',1e2),3),'color',rgb('orangered'),'linewidth',1)
    
    % labels
    axis square tight
    hold off
    xlim(xLims)
    xlabel('time (s)')
    ylabel('BHF LFP (re: baseline)')
    title(deblank(anatomicalLocs{ch2,1}))
    
    
    %% TODO:: plot the stats and active vs passive.
    %% Do I plot both active and passive trials?? - starting with just active.
    subplot(3,2,4)
    title('BART balloon onset stats: ')
    hold on
    maxLim = 24;
    if erH
        text(0,maxLim-2,sprintf('Significant balloon-evoked high gammma increase, p = %0.2f',erP))
    else
        text(0,maxLim-2,sprintf('Insignificant balloon-evoked high gammma increase, p = %0.2f',erP))
    end
    if colorP<=0.05
        text(0,maxLim-4,sprintf('Significant differences among balloon colors \n(Kruskal-Wallis test, p = %0.2f, X^2(%d) = %0.2f)',colorP,Tbl{2,3},Tbl{2,5}))
    else
        text(0,maxLim-4,sprintf('Insignificant differences among balloon colors \n(Kruskal-Wallis test, p = %0.2f, X^2(%d) = %0.2f)',colorP,Tbl{2,3},Tbl{2,5}))
    end
    
    hold off
    
    % deets
    ylim([0 maxLim+1])
    axis off
    
    
    % YELLOW spectra
    ax3 = subplot(3,2,3);
    surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==1),3))),'edgecolor','none');
    caxis([0.5 2])
    set(gca,'Yscale','log');
    view(2)
    axis xy tight square
    colormap(ax3,turbo)
    h = colorbar;
    ylabel(h, 'power re: baseline')
    xlim(xLims)
    ylim([0 params.fpass(2)-50])
    xlabel('time relative to balloon (s)')
    ylabel('frequency (Hz)')
    title(sprintf('n = %d yellow balloons - %s',sum(balloonType==1),deblank(trodeLabels{ch2})))
    
    % ORANGE spectra
    ax4 = subplot(3,2,5);
    surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==2),3))),'edgecolor','none');
    caxis([0.5 2])
    set(gca,'Yscale','log');
    view(2)
    axis xy tight square
    colormap(ax4,turbo)
    h = colorbar;
    ylabel(h, 'power re: baseline')
    xlim(xLims)
    ylim([0 params.fpass(2)-50])
    xlabel('time relative to balloon (s)')
    zlabel('power relative to baseline')
    title(sprintf('n = %d orange balloons - %s',sum(balloonType==2),deblank(trodeLabels{ch2})))
    
    % RED spectra
    ax3 = subplot(3,2,6);
    surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==3),3))),'edgecolor','none');
    caxis([0.5 2])
    set(gca,'Yscale','log');
    view(2)
    axis xy tight square
    colormap(ax3,turbo)
    h = colorbar;
    ylabel(h, 'power re: baseline')
    xlim(xLims)
    ylim([0 params.fpass(2)-50])
    xlabel('time relative to balloon (s)')
    ylabel('frequency (Hz)')
    title(sprintf('n = %d red balloons - %s',sum(balloonType==3),deblank(trodeLabels{ch2})))
    
    % GRAY spectra
    %     ax4 = subplot(3,4,12);
    %     surf(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(abs(mean(Sft(:,(tSec>=-2 & tSec<=2),balloonType==4),3))),'edgecolor','none');
    %     set(gca,'Yscale','log');
    %     view(2)
    %     axis xy tight square
    %     colormap(ax4,parula)
    %     h = colorbar;
    %     ylabel(h, 'power re: baseline')
    %     xlim(xlims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to balloon (s)')
    %     zlabel('power relative to baseline')
    %     title(sprintf('n = %d gray balloons - %s',sum(balloonType==4),deblank(trodeLabels{ch2})))
    
    
    %     % loading phase colormap
    %     load('/home/user1/code/matlab/analyzeMSIT/dependencies/phaseColormap.mat')
    %
    %     % YELLOW phases
    %     ax5 = subplot(3,4,9);
    %     surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(mean(Pft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==1),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==1),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax5,phaseCmap)
    %     xlim(xLims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to balloon (s)')
    %     ylabel('frequency (Hz)')
    %
    %     % ORANGE phases
    %     ax6 = subplot(3,4,10);
    %     surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(mean(Pft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==2),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==2),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax6,phaseCmap)
    %     xlim(xLims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to balloon (s)')
    %
    %     % RED phases
    %     ax5 = subplot(3,4,11);
    %     surf(tSec(tSec>=xLims(1) & tSec<=xLims(2)),scaleFreqs,squeeze(mean(Pft(:,(tSec>=xLims(1) & tSec<=xLims(2)),balloonType==3),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==1),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax5,phaseCmap)
    %     xlim(xLims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to balloon (s)')
    %     ylabel('frequency (Hz)')
    
    % GRAY phases
    %     ax6 = subplot(3,4,16);
    %     surf(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),balloonType==4),3)),'edgecolor','none');
    %     caxis([-pi pi])
    %     set(gca,'Yscale','log');
    %     view(2)c
    %     %     imagesc(tSec(tSec>=-2 & tSec<=2),scaleFreqs,squeeze(mean(Pft(:,(tSec>=-2 & tSec<=2),outcomeType==2),3)),[-pi pi]);
    %     axis xy tight square
    %     h = colorbar;
    %     ylabel(h, 'mean LFP phase (rad)')
    %     colormap(ax6,phaseCmap)
    %     xlim(xlims)
    %     ylim([0 params.fpass(2)-50])
    %     xlabel('time relative to balloon (s)')
    
    
    % overall deets
    halfMaximize(ch2,'left')
    
    % saving data and figures
    saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s/',ptID);
    if ~exist(saveDir,'dir')
        [~,msg] = mkdir(saveDir)
    end
    % print([saveDir '/' sprintf('pt%s_%s_banksvspops_spectrograms_%s',ptID,deblank(trodeLabels{ch2}),normalization)],'-dpdf','-fillpage')
    saveas(ch2,[saveDir '/' sprintf('pt%s_%s_balloonOnset_spectrograms_%s',ptID,deblank(trodeLabels{ch2}),normalization) '.pdf']);
    
    close(ch2)
    
end

% saving staistics
saveDir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
save(fullfile(saveDir,sprintf('%s_balloonOnsetHighGammaStats.mat',ptID)),'balloonOnsetStats')


% do statistics for each channel

% visualize the data here for each channel







