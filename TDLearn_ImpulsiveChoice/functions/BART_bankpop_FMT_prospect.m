function [] = BART_bankpop_FMT_propsect(ptID,whichEEG)
% BART_BANKPOP_FMT_PROSPECT examines reward prediction of FMTheta
%
%   [] = BART_bankpop_FMT_prospect(ptID) analyzes correlations
%   between (frontocentral) EEG and ECoG theta power and BHF activity
%   across trials using linear regression. 
%

%TODO:: just works with pre-cue baseline normalization at this point.

% author: EHS20181005

% ptID = '202004';

nevList = dir(sprintf('~/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% loading behavioral matFile
matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
selectedChans = find(isECoG); % just ECoG, because I'm pulling out EEG manually.
% if strcmp(ptID,'201810'); selectedChans(1:10)  = []; end;
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

% % common average rereferencing w/ PCA [20190905]
% data2K(selectedChans,:) = remove1stPC(data2K(selectedChans,:));

% determineing which EEG electrode to compare against
if nargin<2
whichEEG = 'Fz';
end
whichTrode = contains(trodeLabels,whichEEG);


% downsamplin that trode.
eegChans = find(isEEG);
for eg = 1:length(eegChans)
	if notchFilter
		try
			[b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
			tmp2K(eg,:) = filtfilt(b,a,resample(double(NSX.Data(eg,:)),Fnew,Fs));
		catch
			tmp2K(eg,:) = resample(double(NSX.Data(eg,:)),Fnew,Fs);
		end
	else
		tmp2K(eg,:) = resample(double(NSX.Data(eg,:)),Fnew,Fs);
	end
end
clear NSX NEV
Fs = Fnew;

% removing the common mode across ECoG channels and 
eeg2K_all = remove1stPC(tmp2K);
eeg2K = eeg2K_all(find(isEEG)==find(whichTrode),:);

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
pointsPerTrial = [data.points];


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

	EEGtheta = squeeze(mean(Sft_eeg(scaleFreqsEEG>4 & scaleFreqsEEG<8,:,:)));
	meanBank = mean(EEGtheta(:,outcomeType==1),2);
	errBank = std(EEGtheta(:,outcomeType==1),[],2)./sqrt(sum(outcomeType==1));
	meanPop = mean(EEGtheta(:,outcomeType==2),2);
	errPop = std(EEGtheta(:,outcomeType==2),[],2)./sqrt(sum(outcomeType==2));

	% ECoGtheta = squeeze(mean(mean(Sft(scaleFreqsECoG>4 & scaleFreqsECoG<8,tSec>=0 & tSec<=1,:))));
    % ECoG_BHF = squeeze(mean(mean(Sft(scaleFreqsECoG>70 & scaleFreqsECoG<150,tSec>=0 & tSec<=1,:))));

	% Summary statistics for trials
	nBanks = sum(outcomeType==1);
	nPops = sum(outcomeType==2);
	sortedOutcomeCmap = [repmat(rgb('forestgreen'),sum(outcomeType==1),1); repmat(rgb('orangered'),sum(outcomeType==2),1)];
 
end

% visualize the data here for each channel
if ishandle(ch2); close(ch2); end
figure(ch2)

%% plotting average high gamma responses
subplot(3,2,1)
hold on
patch([tSec fliplr(tSec)],[meanPop'-errPop' fliplr(meanPop'+errPop')],rgb('orangered'),'facealpha',0.5,'edgecolor','none')
plot(tSec,meanPop','color',rgb('orangered'))
text(-2,max(meanPop'+errPop')-1,sprintf('N = %d',uint16(nPops)),'color',rgb('orangered'))
patch([tSec fliplr(tSec)],[meanBank'-errBank' fliplr(meanBank'+errBank')],rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
plot(tSec,meanBank','color',rgb('forestgreen'))
text(-2,max(meanPop'+errPop'),sprintf('N = %d',uint16(nBanks)),'color',rgb('forestgreen'))
%     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
hold off

% plot one deets
axis tight square
xlim([-2 2])
xlabel('time relative to outcome (s)')
ylabel('theta power')
title([deblank(trodeLabels{ch2}) ' -- ' deblank(anatomicalLocs{ch2})])

%% plotting all trials high gamma
subplot(3,2,2)
[~,sortedIdcs] = sort(outcomeType);
hold on
imagesc(tSec,1:nTrials,squeeze(EEGtheta(:,sortedIdcs))')
scatter(repmat(2,1,nTrials),1:nTrials,8,flipud(sortedOutcomeCmap),'filled','s')
hold off

% plot two deets
colorbar
axis tight square xy
xlim([-2 2])
xlabel('time relative to outcome (s)')
ylabel('trials')
title('per trial high gamma')


%% high gamma vs outcome
subplot(3,2,3)
% first fit the data
lineBuffer = 20;

lmPop = fitlm(-pointsPerTrial(outcomeType==2),squeeze(mean(EEGtheta(tSec>0.25 & tSec<1.25,outcomeType==2))),'y ~ 1 + x1');
%     qmPop = fitlm(-pointsPerTrial(outcomeType==2),squeeze(mean(HGmat(ch2,tSec>0.25 & tSec<1.25,outcomeType==2),2))','y ~ 1 + x1^2');
fPopLine = polyval(fliplr(lmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);
%     fPopQuad = polyval(fliplr(qmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);

lmBank = fitlm(pointsPerTrial(outcomeType==1),squeeze(mean(EEGtheta(tSec>0.5 & tSec<2,outcomeType==1))),'y ~ 1 + x1');
%     qmBank = fitlm(pointsPerTrial(outcomeType==1),squeeze(mean(HGmat(ch2,tSec>0.5 & tSec<2,outcomeType==1),2))','y ~ 1 + x1^2');
fBankLine = polyval(fliplr(lmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);
%     fBankQuad = polyval(fliplr(qmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);


% then plot scatters and fitted lines...
hold on
scatter([pointsPerTrial(outcomeType==1) -pointsPerTrial(outcomeType==2)],...
[squeeze(mean(EEGtheta(tSec>0.5 & tSec<2,outcomeType==1))) squeeze(mean(EEGtheta(tSec>0.25 & tSec<1.25,outcomeType==2)))],...
5,sortedOutcomeCmap,'filled')
% plot linear fits
plot(min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer,fPopLine,'color',rgb('orangered'))
plot(min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer,fBankLine,'color',rgb('forestgreen'))

% plot quadratic fits
%     plot(min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer,fPopQuad,'color',rgb('orangered'),'linestyle','--')
%     plot(min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer,fBankQuad,'color',rgb('forestgreen'),'linestyle','--')

hold off

% plot three deets
title('per trial high gamma vs. points gained-lost')
axis square tight

colormap(turbo)
halfMaximize(ch2,'left')


%% statistics
subplot(3,2,4)
hold on
betterBoxplot(1,lmPop.Residuals.Raw.^2,rgb('orangered'),4)
betterBoxplot(2,lmBank.Residuals.Raw.^2,rgb('forestgreen'),4)
hold off

title('Squared residuals from each model')
xlim([0.5 2.5])
axis square

% printing pop stats
subplot(3,2,5)
hold on
text(0,1,'Pop Model:','fontweight','bold')
text(0,0.8,sprintf('intercept: %.2f',lmPop.Coefficients.Estimate(1)))
text(0,0.6,sprintf('slope: %.2f',lmPop.Coefficients.Estimate(2)))
if lmPop.Coefficients.pValue(2)<0.05
	text(0,0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',lmPop.DFE,lmPop.Coefficients.tStat(2),lmPop.Coefficients.pValue(2)))
else
	text(0,0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',lmPop.DFE,lmPop.Coefficients.tStat(2),lmPop.Coefficients.pValue(2)))
end
axis off
hold off

% saving RPE statistics in a structure
prospectStats(ch2).trodeLabel = deblank(trodeLabels{ch2});
prospectStats(ch2).PopLinearModel = lmPop;
prospectStats(ch2).BankLinearModel = lmBank;

% saving figures
saveDir = sprintf('~/Dropbox',ptID);
if exist(fullfile(saveDir,'FMTprospect'),'dir')
	saveas(ch2,fullfile(saveDir,'FMTprospect',sprintf('pt%s_banksvspops_FMT_prospect.pdf',ptID)))
else
	mkdir(saveDir,'FMTprospect')
	saveas(ch2,fullfile(saveDir,'FMTprospect',sprintf('pt%s_banksvspops_FMT_prospect.pdf',ptID)))
end

close(ch2)

end



