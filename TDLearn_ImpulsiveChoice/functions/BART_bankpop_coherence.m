function [cohdata] = BART_bankpop_coherence(ptID)
% BART_BANKPOP_COHERENCE analyzes and visualizes pairwise coherence
%
%   [BARTstats] = BART_bankpop_coherence(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   THIS FUNCTION USES WAVELET COHERENCE. Use BART_bankpop_MTcoherence for
%   multitaper.

% author: EHS20181005

% in order to run a single patient...
ptID = '201905';


%% loading data
nevList = dir(sprintf('~/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end

% trode details.
[trodeLabels,isECoG,isEEG,~,anatomicalLocs] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% channel specification; removing whites.
tmpChans = find(isECoG);
whites = contains(anatomicalLocs,'White');
if strcmp(ptID,'201810'); whites(1:10) = true; end
selectedChans = tmpChans(~whites);
nChans = length(selectedChans);
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;


%% resampling LFP at Fnew sampling frequency
notchFilter = true;
Fnew = 250;
for ch = 1:nChans
    if notchFilter
        [b,a] = iirnotch(60/(Fnew/2),(60/(Fnew/2))/50);
        data2K(ch,:) = filtfilt(b,a,resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs));
    else
        data2K(ch,:) = resample(double(NSX.Data(selectedChans(ch),:)),Fnew,Fs);
    end
end
clear NSX NEV
Fs = Fnew;

% % common average re-referencing with 1st PC.
% data2K = remove1stPC(data2K);

% timing parameters.
pre = 3;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpass = [1 ceil(Fnew./2)]; % pick a value that's 50 higher for wavelets
params.pad = 0;
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
            if (ppData.Event.trigTimes(1)./3e4)<secsPreTask
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

% Summary statistics for each channel
nBanks = sum(outcomeType==1);
nPops = sum(outcomeType==2);


myAnalyses = true;
if myAnalyses
    % epoching data
	fprintf('\nepoching data...\n');
    LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
    % plotting bank/pop responses
    for ch2 = 1:nChans
        % epoch the spectral data for each channel.\
        for tt = 1:nTrials
            % epoch the data here [channels X samples X trials]
            LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
            
            % do spectral claculations here
            [tmpW(:,:,tt,ch2),period,scale] = basewaveERP(LFPmat(ch2,:,tt),Fs,params.fpass(1),params.fpass(2),6,0);
        end
    end
    clear data2K
    
    % determining spectral scales in Hz
    tmpF = period.^-1;
    
    % trimming W to make these analyses possible
    W = tmpW(tmpF<params.fpass(2)-50,tSec>-pre+1 & tSec<post-1,:,:);
    tSec = tSec(tSec>-pre+1 & tSec<post-1);
    fHz = tmpF(tmpF<params.fpass(2)-50);
    clear tmpW tmpF
    
    % pairs
    cohdata = struct();
    pairs = nchoosek(selectedChans,2);
    cohdata.pairs = pairs;
    for pwchs = 1:length(pairs)
        updateUser('coherence calculation for pair',pwchs,100,length(pairs));
        % doing  coherence calculation       
        for trl = 1:length(outcomeType)
	        XS = W(:,:,trl,pairs(pwchs,1)).*conj(W(:,:,trl,pairs(pwchs,2)));
		    tmpC = ((XS).^2)./(((W(:,:,trl,pairs(pwchs,1))).^2).* ((W(:,:,trl,pairs(pwchs,2)))).^2);
 			C(:,:,trl) = abs(tmpC);
			% C(:,:,trl) = imgaussfilt(abs(real(squeeze(tmpC))),2);
        end
        clear tmpC
        
        %% permutation stats
        % true value
	    Csize = length(tSec)*length(fHz);
        for ts=1:Csize
            [i(ts),j(ts)] = ind2sub([size(C,1),size(C,2)],ts);
            [pData(ts),~,statsData{ts}] = ranksum(squeeze(C(i(ts),j(ts),outcomeType==1)),squeeze(C(i(ts),j(ts),outcomeType==2)));
        end
        
        % permutation distribution
        cohData.nPerms = 1;
        for prm = 1:cohData.nPerms
			updateUser('did permutation stats',prm,100,cohData.nPerms);
            % make random grouping vector
            groupingVec(prm,:) = outcomeType(randperm(length(outcomeType)))==1;
            
            % run test
            for ts=1:Csize
                [i,j] = ind2sub([size(C,1),size(C,2)],ts);
                [cohData.pPerm(ts,prm),~,cohData.statsPerm{ts}] = ranksum(squeeze(C(i,j,groupingVec(prm,:))),squeeze(C(i,j,~groupingVec(prm,:))));
            end
            
            if prm==cohData.nPerms
				% final variables added to cohdata
				cohData.trialPermutations = groupingVec;
				cohData.C = C;

                % save confidence image for comparison. TODO::: Can I save it to the server?
				saveDir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
				fName = sprintf('%s_FFcoherence_ch%dand%d.mat',ptID,pairs(pwchs,1),pairs(pwchs,2));
				save(fullfile(saveDir,fName),'cohData','-v7.3')
            end
        end
        
        % plotting for each channel pair.
        figure(pwchs)
        
        % first graph: channel 1 banks
        subplot(3,2,1)
        surf(tSec,fHz,mean((abs(W(:,:,outcomeType==1,pairs(pwchs,1))))./repmat(mean(mean((abs(W(:,tSec>-2 & tSec<-1,:,pairs(pwchs,1)))),3),2),[1 length(tSec) nBanks]),3),'edgecolor','none')
        set(gca,'Yscale','log');
        view(2)
        axis xy tight square
        xlim([-pre+2 post-1])
        zlim([0.5 2])
        colorbar
        title(sprintf('bank spectrum: %s',anatomicalLocs{pairs(pwchs,1)}))
        
        % second graph: channel 1 pops
        subplot(3,2,2)
        surf(tSec,fHz,mean((abs(W(:,:,outcomeType==2,pairs(pwchs,2))))./repmat(mean(mean((abs(W(:,tSec>-2 & tSec<-1,:,pairs(pwchs,2)))),3),2),[1 length(tSec) nPops]),3),'edgecolor','none')
        set(gca,'Yscale','log');
        view(2)
        axis xy tight square
        xlim([-pre+2 post-1])
        zlim([0.5 2])
        colorbar
        title(sprintf('pop spectrum: %s',anatomicalLocs{pairs(pwchs,2)}))
        
        % third graph: channel 2 banks
        subplot(3,2,3)
        surf(tSec,fHz,mean((abs(W(:,:,outcomeType==1,pairs(pwchs,1))))./repmat(mean(mean((abs(W(:,tSec>-2 & tSec<-1,:,pairs(pwchs,1)))),3),2),[1 length(tSec) nBanks]),3),'edgecolor','none')
        set(gca,'Yscale','log');
        view(2)
        axis xy tight square
        xlim([-pre+2 post-1])
        colorbar
        title(sprintf('bank spectrum: %s',anatomicalLocs{pairs(pwchs,1)}))
        
        % fourth graph: channel 1 pops
        subplot(3,2,4)
        surf(tSec,fHz,mean((abs(W(:,:,outcomeType==2,pairs(pwchs,2))))./repmat(mean(mean((abs(W(:,tSec>-2 & tSec<-1,:,pairs(pwchs,2)))),3),2),[1 length(tSec) nPops]),3),'edgecolor','none')
        set(gca,'Yscale','log');
        view(2)
        axis xy tight square
        xlim([-pre+2 post-1])
        zlim([0.5 2])
        colorbar
        title(sprintf('pop spectrum: %s',anatomicalLocs{pairs(pwchs,2)}))
        
        % 5th graph: bank coherence
        subplot(3,2,5)
        surf(tSec,fHz,squeeze(mean(C(:,:,outcomeType==1),3)),'edgecolor','none')
        set(gca,'Yscale','log');
        axis xy tight square
        xlim([-pre+2 post-1])
        zlim([0 0.5])
        view(2)
        colorbar
        title('coherence between above channels - banks')
        
        % 6th graph: pop coherence
        subplot(3,2,6)
        surf(tSec,fHz,squeeze(mean(C(:,:,outcomeType==2),3)),'edgecolor','none')
        set(gca,'Yscale','log');
        axis xy tight square
        xlim([-pre+2 post-1])
        zlim([0 0.5])
        view(2)
        colorbar
        title('coherence between above channels - pops')
        
        % turbo colormap
        colormap(turbo)
        
        halfMaximize(pwchs,'left')
        saveas(pwchs,fullfile(saveDir,sprintf('%s_FieldFieldCoherenceChans.pdf',ptID,ptID)))
        
        keyboard

		close(pwchs)
    end
    
    % save(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/%s_pairwiseCoherence.mat',ptID,ptID),'cohdata','-v7.3')
    
end







% these are wavelet coherograms based on the wavelet coherence toolboc from grinsted
geoAnalyses = false;
if geoAnalyses
    
    %% epoching data
    LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
    % plotting bank/pop responses
    for ch2 = 1:nChans
        % epoch the spectral data for each channel.\
        for tt = 1:nTrials
            % epoch the data here [channels X samples X trials]
            LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
            
        end
    end
    
    %% figuring out which pairs to check out and setting up labels
    pairs = nchoosek(selectedChans,2);
    for pwchs = 1:length(pairs)
        seriesname{1} = anatomicalLocs{pairs(pwchs,1)};
        seriesname{2} = anatomicalLocs{pairs(pwchs,2)};
        
        
        
        %% Continuous wavelet transform (CWT)
        % The CWT expands the time series into time
        % frequency space.
        figure('color',[1 1 1])
        tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];
        subplot(2,2,1);
        wt(d1);
        title(seriesname{1});
        set(gca,'xlim',tlim);
        subplot(2,2,2)
        wt(d2)
        title(seriesname{2})
        set(gca,'xlim',tlim)
        
        
        %% Cross wavelet transform (XWT)
        % The XWT finds regions in time frequency space where
        % the time series show high common power.
        subplot(2,2,3)
        xwt(d1,d2)
        title(['XWT: ' seriesname{1} '-' seriesname{2} ] )
        
        
        %% Wavelet coherence (WTC)
        % The WTC finds regions in time frequency space where the two
        % time series co-vary (but does not necessarily have high power).
        subplot(2,2,4)
        wtc(d1,d2)
        title(['WTC: ' seriesname{1} '-' seriesname{2} ] )
        
    end
end
