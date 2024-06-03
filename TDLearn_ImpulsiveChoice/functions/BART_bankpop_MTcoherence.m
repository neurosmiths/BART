function [cohdata] = NBACK_Hit_MTcoherence(ptID)
% NBACK_HIT_COHERENCE analyzes and visualizes pairwise coherence
%
%   [NBACKstats] = NBACK_Hit_coherence(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20181005

% in order to run a single patient...
 ptID = '202202';

nevList = dir(['D:\Data\preProcessed\NBACK_preprocessed\' ptID '\Data\' ptID '_N_Back_*.txt']);
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
[trodeLabels,isECoG,isEEG,~,anatomicalLocs] = ptTrodesNBACK(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns6']));

% channel specification; removing whites. 
tmpChans = find(isECoG);
whites = contains(anatomicalLocs,'White');
if strcmp(ptID,'202201'); whites(1:10) = true; end
selectedChans = tmpChans(~whites);
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
clear NSX NEV
Fs = Fnew;

% % common average re-referencing with 1st PC.
% data2K = remove1stPC(data2K);

% timing parameters.
pre = 2;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% TF parameters
params.fpass = [1 90]; % pick a value that's 50 higher for wavelets
params.tapers = [7 13];
params.pad = 0;
params.Fs = Fs;
params.dBconversion = false;
params.normalized = false;
params.theoreticalNorm = false;
params.baseline = true;
params.err = [2 0.05];
params.trialave = 1;
movingwin = [0.8 0.05];

% picking pre-bird baseline period
if params.baseline
    baselineType = 'preCue';
    switch baselineType
        case {'preCue'}
            bP = [-2.5 -1.5];
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
stimTimes = trigTimes(trigs==2);
respTimes = trigTimes(trigs==3 | trigs==4);

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% task parameters in chronological order..
outcomeType = trigs(sort([find(trigs==3); find(trigs==4)])); % 3 = hits, 4 = misses
nTrials = length(outcomeType);

% Summary statistics for each channel
nHits = sum(outcomeType==3);
nMiss = sum(outcomeType==4);

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
% plotting hit/miss responses
for ch2 = 1:nChans    
    % epoch the spectral data for each channel.\
    for tt = 1:nTrials
        % epoch the data here [channels X samples X trials]
        LFPmat(ch2,:,tt) = data2K(ch2,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    end    
end


% pairs
cohdata = struct();
pairs = nchoosek(1:length(selectedChans),2);
cohdata.pairs = pairs;
for pwchs = 1:length(pairs)    
    % doing  coherence calculation
    tic;
    fprintf('\ncalculating coherence with jacknife error... \n(%d out of %d channel pairs)...',pwchs,length(pairs))
    [cohdata(pwchs).Cbank,~,cohdata(pwchs).S12bank,cohdata(pwchs).S1bank,cohdata(pwchs).S2bank,cohdata(pwchs).t,cohdata(pwchs).f,cohdata(pwchs).confCbank,cohdata(pwchs).phiSTDbank,cohdata(pwchs).Cerrbank]...
        = cohgramc(squeeze(LFPmat(pairs(pwchs,1),:,outcomeType==1)),squeeze(LFPmat(pairs(pwchs,2),:,outcomeType==1)),movingwin,params);
    [cohdata(pwchs).Cpop,~,cohdata(pwchs).S12pop,cohdata(pwchs).S1pop,cohdata(pwchs).S2pop,cohdata(pwchs).t,cohdata(pwchs).f,cohdata(pwchs).confCpop,cohdata(pwchs).phiSTDpop,cohdata(pwchs).Cerrpop]...
        = cohgramc(squeeze(LFPmat(pairs(pwchs,1),:,outcomeType==2)),squeeze(LFPmat(pairs(pwchs,2),:,outcomeType==2)),movingwin,params);
    A = toc;
    fprintf('\n...finished. Took %.2f minutes.',A/60)
    t = cohdata(pwchs).t-pre;
    
    % plotting
    tic;
    figure(pwchs)
    
    subplot(3,2,1)
    imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).S1bank./repmat(mean(cohdata(pwchs).S1bank(t<-1,:)),length(t),1))')
    axis xy
    colorbar
    title(sprintf('baseline norm spectrum [banks]: %s',anatomicalLocs{pairs(pwchs,1)}))
    
    subplot(3,2,2)
    imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).S1pop./repmat(mean(cohdata(pwchs).S1pop(t<-1,:)),length(t),1))')
    axis xy
    colorbar
    title(sprintf('baseline norm spectrum [pops]: %s',anatomicalLocs{pairs(pwchs,2)}))

    subplot(3,2,3)
    imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).S2bank./repmat(mean(cohdata(pwchs).S2bank(t<-1,:)),length(t),1))')
    axis xy
    colorbar
    title(sprintf('baseline norm spectrum [banks]: %s',anatomicalLocs{pairs(pwchs,1)}))
    
    subplot(3,2,4)
    imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).S2pop./repmat(mean(cohdata(pwchs).S2pop(t<-1,:)),length(t),1))')
    axis xy
    colorbar
    title(sprintf('baseline norm spectrum [pops]: %s',anatomicalLocs{pairs(pwchs,2)}))

    subplot(3,2,5)
%     imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).C./repmat(mean(cohdata(pwchs).C(t<-1,:)),length(t),1))')
    imagesc(t,cohdata(pwchs).f,cohdata(pwchs).Cbank')
    colorbar
    axis xy
    title('coherence between two above channels')

    subplot(3,2,6)
%     imagesc(t,cohdata(pwchs).f,(cohdata(pwchs).C./repmat(mean(cohdata(pwchs).C(t<-1,:)),length(t),1))')
    imagesc(t,cohdata(pwchs).f,cohdata(pwchs).Cpop')
    colorbar
    axis xy
    title('coherence between two above channels')

    colormap(turbo)
    
    halfMaximize(pwchs,'left')
        
    saveas(pwchs,sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/%s_MTFieldFieldCoherenceChans_%dand%d.pdf',ptID,ptID,pairs(pwchs,1),pairs(pwchs,2)))
    A = toc; 
    fprintf('\n...finished plotting and saving. Took %.2f minutes.',A/60)

    close(pwchs)

    save(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/%s_pairwiseCoherence.mat',ptID,ptID),'cohdata','-v7.3')
end


