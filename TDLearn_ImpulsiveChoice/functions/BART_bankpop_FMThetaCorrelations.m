function [] = BART_bankpop_FMThetaCorrelations(ptID,whichEEG)
% BART_BANKPOP_FMTHETACORRELATIONS analyzes correlations betewen EEG & ECoG
%
%   [] = BART_bankpop_FMThetaCorrelations(ptID) analyzes correlations
%   between (frontocentral) EEG and ECoG theta power and BHF activity
%   across trials using linear regression. 
%

%TODO:: just works with pre-cue baseline normalization at this point.

% author: EHS20181005

 ptID = '202004';
 
parentDir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev'];
nevList = dir(parentDir);
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
selectedChans = find(isECoG); % just ECoG, because I'm pulling out EEG manually.
nChans = length(selectedChans);
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% resampling LFP at Fnew sampling frequency
notchFilter = true;
Fnew = 500;
for ch = 1:nChans
    if notchFilter
        [b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        data2K(ch,:) = filtfilt(b,a,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end

% common average rereferencing w/ PCA [20190905]
% data2K(selectedChans,:) = remove1stPC(data2K(selectedChans,:));

% determineing which EEG electrode to compare against
if nargin<2
whichEEG = 'Fz';
end
whichTrode = contains(trodeLabels,whichEEG);

% downsamplin that trode.
if notchFilter
    try
        [b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        eeg2K = filtfilt(b,a,resample(double(NSX.Data(whichTrode,:)),Fnew,Fs));
    catch
        eeg2K = resample(double(NSX.Data(whichTrode,:)),Fnew,Fs);
    end
else
    eeg2K = resample(double(NSX.Data(whichTrode,:)),Fnew,Fs);
end
clear NSX NEV
Fs = Fnew;

% timing parameters.
pre = 3;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpassECoG = [1 200]; % pick a value that's ~50 higher for wavelets
params.fpassEEG = [1 100]; % pick a value that's ~50 higher for wavelets
params.Fs = Fs;
% normalization parameters
params.dBconversion = false;
params.normalized = false;
params.theoreticalNorm = false;
params.baseline = true;

% picking pre-stimulus baseline period
if params.baseline
    baselineType = 'preCue';
    switch baselineType
        case {'preCue'}
            bP = [-1.5 -0.5];
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
respTimes = trigTimes(trigs==23);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);
trialMap = zeros(nTrials,3);
trialMap(outcomeType==1,:) = repmat(rgb('forestgreen'),sum(outcomeType==1),1); 
trialMap(outcomeType==2,:) = repmat(rgb('orangered'),sum(outcomeType==2),1);  


% epoching  data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting bank/pop responses
for ch2 = 1:nChans
    % epoch the spectral data for each channel.\
    for tt = 1:nTrials
        fprintf('\ndoing spectral calculations for trial %d and channel %d...',tt,ch2)
        % epoch ECoG here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*respTimes(tt))-Fs*pre:floor(Fs*respTimes(tt))+Fs*post);
        
        % epoching EEG data
        EEGmat(:,tt) = eeg2K(floor(Fs*respTimes(tt))-Fs*pre:floor(Fs*respTimes(tt))+Fs*post);
        
        % do spectral claculations here
        [W,periodECoG,scale] = basewaveERP(LFPmat(ch2,:,tt),Fs,params.fpassECoG(1),params.fpassECoG(2),6,0);
        [Weeg,periodEEG,~] = basewaveERP(EEGmat(:,tt),Fs,params.fpassEEG(1),params.fpassEEG(2),6,0);
        
        if params.dBconversion
            normalization = 'decibels';
            Sft(:,:,tt) = abs(10*log10((W)));
        elseif params.baseline
            Sft(:,:,tt) = abs(W);
            Sft_eeg(:,:,tt) = abs(Weeg);
        elseif params.normalized
            normalization = 'normByFreq';
            Sft(:,:,tt) = normlogspec(abs((W))')';
        elseif params.theoreticalNorm
            normalization = 'theoreticalNorm';
            Sft(:,:,tt) = abs((W))...
                ./repmat(1./linspace(params.fpassECoG(1),params.fpassECoG(2),length(scale))',1,size(W,2));
        end
        
        % phase angle.
        Pft(:,:,tt) = angle(W);
        
    end
    
    % doing baseline normalization
    if (params.baseline & baselineType=='preCue')
        normalization = [baselineType 'baselineNorm'];
        tmp = Sft./repmat(mean(mean(Sft(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft,2),size(Sft,3));
        Sft = tmp;
        Sft_eeg = Sft_eeg./repmat(mean(mean(Sft_eeg(:,tSec>bP(1) & tSec<bP(2),:),2),3),1,size(Sft_eeg,2),size(Sft_eeg,3));
        clear tmp
    elseif (params.baseline & baselineType=='preTask')
        normalization = [baselineType 'baselineNorm'];
        % TODO:: add pre task baseline normalization.
    end
    
    % determining spectral scales in Hz
    scaleFreqsECoG = periodECoG.^-1;
    scaleFreqsEEG = periodEEG.^-1;
    
    % Summary statistics for trials
    nBanks = sum(outcomeType==1);
    nPops = sum(outcomeType==2);
    
    % TODO:: run statistical models here
    % simple linear model relating theta from a ECoG contact to theta from
    % an EEG contact.
    EEGtheta = squeeze(mean(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,tSec>=0 & tSec<=1,:))));
    ECoGtheta = squeeze(mean(mean(Sft(scaleFreqsECoG>4 & scaleFreqsECoG<8,tSec>=0 & tSec<=1,:))));
    ECoG_BHF = squeeze(mean(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=0 & tSec<=1,:))));
    
    % linear model
    trodeLMtheta = fitlm(EEGtheta,ECoGtheta,'linear','intercept',true);
    trodeLM_BHF = fitlm(EEGtheta,ECoG_BHF,'linear','intercept',true);
    
    % TODO:: save statistics
    % saving data and figures
    EEGECoGtrodeLM.ptID = ptID;
    EEGECoGtrodeLM.EEGtrode = whichEEG;
    EEGECoGtrodeLM.ECoGtrode = deblank(trodeLabels{selectedChans(ch2)});
    EEGECoGtrodeLM.trodeLMtheta = trodeLMtheta;
    EEGECoGtrodeLM.trodeLM_BHF = trodeLM_BHF;
    dataDir = sprintf('D:/Data/Rhiannon/BART_RLDM_outputs/BatchAnalyze/FMT',ptID);
    save([dataDir '/' sprintf('pt%s_%s_outcome_%s_EEGECoGcorrelations_%s.mat'...
        ,ptID,whichEEG,deblank(trodeLabels{selectedChans(ch2)}),normalization)],'EEGECoGtrodeLM')
    
    SIG = true;
    if SIG
        % visualize the data here for EEG
        if ishandle(ch2); close(ch2); end
        figure(ch2)
        colormap(parula)
        xlims = [-1 2];
        
%         % pop spectra EEG
%         ax4 = subplot(3,2,1);
%         surf(tSec(tSec>=-2 & tSec<=2),scaleFreqsEEG,squeeze(abs(mean(Sft_eeg(:,tSec>=-2 & tSec<=2,outcomeType==2),3))),'edgecolor','none');
%         set(gca,'Yscale','log');
%         view(2)
%         axis xy tight square
%         colormap(ax4,parula)
%         h = colorbar;
%         ylabel(h, 'power re: baseline')
%         xlim(xlims)
%         ylim([0 params.fpassEEG(2)-50])
%         xlabel('time relative to outcome (s)')
%         zlabel('power relative to baseline')
%         title([whichEEG ' [balloon pops]'])
%         
%         % pop spectra ECOG
%         ax4 = subplot(3,2,2);
%         surf(tSec(tSec>=-2 & tSec<=2),scaleFreqsECoG,squeeze(abs(mean(Sft(:,tSec>=-2 & tSec<=2,outcomeType==2),3))),'edgecolor','none');
%         set(gca,'Yscale','log');
%         view(2)
%         axis xy tight square
%         colormap(ax4,parula)
%         h = colorbar;
%         ylabel(h, 'power re: baseline')
%         xlim(xlims)
%         ylim([0 params.fpassECoG(2)-50])
%         xlabel('time relative to outcome (s)')
%         zlabel('power relative to baseline')
%         title([trodeLabels(selectedChans(ch2)) ' [balloon pops]'])
        
%         % theta power from EEG
%         subplot(3,2,3)
%         hold on
%         plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,tSec>=-2 & tSec<=2,outcomeType==1),1),2,'movmean',2e2)),'color',rgb('palegreen'))
%         plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,tSec>=-2 & tSec<=2,outcomeType==2),1),2,'movmean',2e2)),'color',rgb('peachpuff'))
%         plot(tSec(tSec>=-2 & tSec<=2),mean(smoothdata(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,tSec>=-2 & tSec<=2,outcomeType==1),1),2,'movmean',1e2),3),'color',rgb('forestgreen'),'linewidth',1)
%         plot(tSec(tSec>=-2 & tSec<=2),mean(smoothdata(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,tSec>=-2 & tSec<=2,outcomeType==2),1),2,'movmean',1e2),3),'color',rgb('orangered'),'linewidth',1)
%         % labels
%         axis square tight
%         xl = xlim;
%         yl = ylim;
%         text(xl(2)-1,yl(2)-0.2,'banked','color',rgb('forestgreen'),'fontweight','bold')
%         text(xl(2)-1,yl(2)-0.4,'popped','color',rgb('orangered'),'fontweight','bold')
%         hold off
%         xlim(xlims)
%         xlabel('time (s)')
%         ylabel('theta power (re: baseline)')
%         
%         % BHF power from ECoG
%         subplot(3,2,4)
%         hold on
%         plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=-2 & tSec<=2,outcomeType==1),1),2,'movmean',2e2)),'color',rgb('palegreen'))
%         plot(tSec(tSec>=-2 & tSec<=2),squeeze(smoothdata(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=-2 & tSec<=2,outcomeType==2),1),2,'movmean',2e2)),'color',rgb('peachpuff'))
%         plot(tSec(tSec>=-2 & tSec<=2),mean(smoothdata(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=-2 & tSec<=2,outcomeType==1),1),2,'movmean',1e2),3),'color',rgb('forestgreen'),'linewidth',1)
%         plot(tSec(tSec>=-2 & tSec<=2),mean(smoothdata(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=-2 & tSec<=2,outcomeType==2),1),2,'movmean',1e2),3),'color',rgb('orangered'),'linewidth',1)
%         % labels
%         axis square tight
%         xl = xlim;
%         yl = ylim;
%         text(xl(2)-1,yl(2)-0.2,'banked','color',rgb('forestgreen'),'fontweight','bold')
%         text(xl(2)-1,yl(2)-0.4,'popped','color',rgb('orangered'),'fontweight','bold')
%         hold off
%         xlim(xlims)
%         xlabel('time (s)')
%         ylabel('BHF power (re: baseline)')
%         title(trodeLabels(selectedChans(ch2)))
        
        % plot significant regression models here.
        subplot(3,2,5)
        hold on
        scatter(trodeLMtheta.Variables.x1,trodeLMtheta.Variables.y,5,trialMap,'filled')
        thetaH = trodeLMtheta.plot;
        thetaH(1).Color = 'none';
        thetaH(2).Color = [0 0 0];
        thetaH(3).Color = [0 0 0];
        thetaH(4).Color = [0 0 0];
        hold off
        axis tight square
        xlabel('frontal midline theta')
        ylabel('intracranial theta')
        title(sprintf('FDR-c p-value = %.3f; R squared = %.2f',trodeLMtheta.Coefficients.pValue(2)*nChans,trodeLMtheta.Rsquared.Ordinary))
        
        subplot(3,2,6)
        hold on
        scatter(trodeLM_BHF.Variables.x1,trodeLM_BHF.Variables.y,5,trialMap,'filled')
        bhfH = trodeLM_BHF.plot;
        bhfH(1).Color = 'none';
        bhfH(2).Color = [0 0 0];
        bhfH(3).Color = [0 0 0];
        bhfH(4).Color = [0 0 0];
        hold off
        axis tight square
        xlabel('frontal midline theta')
        ylabel('intracranial BHF')
        title(sprintf('FDR-c p-value = %.3f; R squared = %.2f',trodeLM_BHF.Coefficients.pValue(2)*nChans,trodeLM_BHF.Rsquared.Ordinary))
        
        % overall deets
        halfMaximize(ch2,'left')
        
        % saving data and figures
        saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s',ptID);
        if ~exist(saveDir,'dir')
            [~,msg] = mkdir(saveDir)
        end
        print([saveDir '/' sprintf('pt%s_%s_outcome_%s_EEGECoGcorrelations_%s',ptID,whichEEG,deblank(trodeLabels{selectedChans(ch2)}),normalization)],'-dpdf')
        
        close(ch2)
    end
end



