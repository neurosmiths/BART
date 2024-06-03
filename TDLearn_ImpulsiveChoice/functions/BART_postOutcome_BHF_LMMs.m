function [trodeLMM] = BART_balloonOnset_postOutcome_BHF(ptID,BARTdir)
% BART_BALOONONSET_POSTOUTCOME_BHF examine post outcome encoding of risk/reward
%
%   [BARTstats] = BART_balloonOnset_postOutcome_bhf(ptID,nevFile) analyzes LFP data
%   for the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20200823

% ptID = '201901'; % this line only for debugging

if nargin<2
    BARTdir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev'];
end
nevList = dir([BARTdir '\*.nev']);

% loading PTB BHV
load(fullfile(BARTdir,[ptID '.bartBHV.mat']))

% electrode information.
[trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

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
Fs = NSX.MetaTags.SamplingFreq;

% timing parameters.
pre = 2;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

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
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% only including complete trials.
nTrials = length(respTimes);
balloonType = balloonType(1:nTrials);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
isCTRL = isCTRL(1:nTrials);
outcome = ~logical(trigs(trigs==25 | trigs==26)-25); % 1 = bank; 0 = pop
postOutcome = [NaN; outcome(1:end-1)]; % 

% colormap for balloon color. #uselelss???
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% outcome color map
oMap(1,:) = rgb('forestgreen');
oMap(2,:) = rgb('orangered');

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% high notch for one patient with a high noise harmonic
if isequal(ptID,'201903')
    [b120,a120] = iirnotch(120/(Fs/2),60/(Fs/2)/35);
end

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

% aligning on balloon appearances to start.
% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    CUEmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*balloonTimes(tt))-Fs*pre:floor(Fs*balloonTimes(tt))+Fs*post);
    OUTmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        % matrix of high gamma signals. [channels X samples X trtials]
        cueHG(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(CUEmat(ch,:,tt)))));
        outHG(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
end

%% setting up regressors for linear models.
% [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]

% permutation tests between conditions for each channel
if ~exist(sprintf('%s/postOutcomeBHFclusterStats_%s.mat',nevList.folder,ptID),'file')
    tic
    [pval,t_orig,clust_info,~,est_alpha] = clust_perm2(HGmat(:,:,postOutcome==0),HGmat(:,:,postOutcome==1),adjacentChanMat);
    A = toc;
    fprintf('\ncluster-based permutation stats took %d seconds',A)
    save(sprintf('%s/postOutcomeBHFclusterStats_%s.mat',nevList.folder,ptID),'pval','t_orig','clust_info','est_alpha')
else
    fprintf('stats have already been done. Aren"t you lucky!')
    load(sprintf('%s/postOutcomeBHFclusterStats_%s.mat',nevList.folder,ptID),'pval','t_orig','clust_info','est_alpha')
end

% how much do we smooth the data?
% Fs = 1000 for these experiments.
% so this would smooth by 200 ms moving average.
smoothFactor = Fs./5;
% plotting responses
for ch2 = 1:nChans
    % setting up variables for statistics
    statTimeWin = [0 1]; % in seconds.

    % plotting
    if ishandle(ch2); close(ch2); end
    figure(ch2)

    % plotting the risk and reward categories
    subplot(3,2,2)
    hold on
    for rsk = 1:2
	    patch([tSec fliplr(tSec)],[(smooth(mean(HGmat(ch2,:,postOutcome==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,postOutcome==(rsk-1)),[],3),smoothFactor)'./nTrials))...
		    fliplr(smooth(mean(HGmat(ch2,:,postOutcome==(rsk-1)),3),smoothFactor)'-(smooth(std(HGmat(ch2,:,postOutcome==(rsk-1)),[],3),smoothFactor)'./nTrials))]...
			,oMap(rsk,:),'facealpha',0.3,'edgecolor','none')
		plot(tSec,smooth(mean(HGmat(ch2,:,postOutcome==(rsk-1)),3),smoothFactor)','color',oMap(rsk,:))
	end
	plot(tSec,~(clust_info.pos_clust_ids(ch2,:)==0),'r*')
	hold off

	% deets
	hold off
	axis tight square
	xlim([-pre+1 post-1])
	xlabel('time relative to balloon onset (s)')
	ylabel('broadband high frequency power')
	keyboard
	title(sprintf('%s (%s) - red stars: significant clusters -',deblank(trodeLabels{ch2}),char(anatomicalLocs{ch2})));


% 	subplot(3,2,3)
% 	notBoxPlot(riskVar,riskGroups+1);
% 	% plotting pairwise comparisons for risk.
% 	if BARTstats(ch2).pRisk<0.05
% 		pairwiseComps = multcompare(BARTstats(ch2).statsRisk,'Display','off');
% 		sigComps = pairwiseComps(pairwiseComps(:,6)<0.05,:);
% 		hold on
% 		for ls = 1:sum(pairwiseComps(:,6)<0.05)
% 			line([sigComps(ls,1) sigComps(ls,1)],[max(riskVar)+0.4 max(riskVar)+0.4],'linewidth',2)
% 		end
% 	end
% 	axis square
% 	xlabel('increasing risk: none, yellow, orange, red')
% 	ylabel('mean BHF power')
% 	title('lines: significant pairwise comparisons')
% 

	% saving figures
	saveDir = ['D:\Data\Rhiannon\BART_RLDM_output\' ptID ''] %does not exist
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    halfMaximize(ch2,'left')
	saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_postOutcome_balloonOnset_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))

end










