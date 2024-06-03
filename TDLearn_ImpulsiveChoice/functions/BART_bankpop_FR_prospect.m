function [] = BART_bankpop_FR_prospect(ptID)   % BARTstats
% BART_BANKPOP_FIRINGRATES analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_bankpop_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.

% author: EHS20181005
ptID = '202107';

nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end
% [trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% loading behavioral matFile
matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)
pointsPerTrial = [data.points];

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% timing parameters: how muhc data to epoch, and then which window to test. 
pre = 3;
post = 3;
pWin = [0.25 1.75];

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
[~,sortedTrialIdcs] = sort(outcomeType);
nTrials = length(outcomeType);

% overall numbers of outcome types
nBanks = sum(outcomeType==1);
nPops = sum(outcomeType==2);

% standard 3-D [chan, unit, timestamp (seconds)] matrix. 
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets. 
inclChans = unique(ChanUnitTimestamp(:,1));
%inclChans = 113:128;
inclChans = 97:112;
[microLabels,microPts] = microLabelsBART(ptID);
% TODO:: THIS LINE IS BROKEN. IT REQUIRES THE MAGIC NUMBER TO BE THE STARING CHANNEL, SO DOESN'T WORK IF UNITS WERE RECORDED ON CHS 17-32 OF BANK D
% inclChans(inclChans-112>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.
nChans = length(inclChans);

% looping over Channels
for ch = 1:nChans
    % looping over number of units in the AP data
    nUnits = length(unique(ChanUnitTimestamp(inclChans(ch).*ones(size(ChanUnitTimestamp,1),1)==ChanUnitTimestamp(:,1),2)));

	for un = 1:nUnits
		fprintf('\nprocessing and plotting for unit %d of %d',un,nUnits)

        % getting unit times for the current channel and unit.
        unitTimes = ChanUnitTimestamp(ChanUnitTimestamp(:,1)==inclChans(ch) & ChanUnitTimestamp(:,2)==un,3); % in seconds
        
		% visualize the data here for each unit
        figure(un)
        subplot(3,2,1)
        
        % loooping over trials
        for tt = 1:nTrials
            
            % putting the data in a structure
            spikes.channel(ch).unit(un).trial(tt).times = unitTimes(unitTimes>outcomeTimes(tt)-pre & unitTimes<outcomeTimes(tt)+post) - repmat(outcomeTimes(tt)-pre,length(unitTimes(unitTimes>outcomeTimes(tt)-pre & unitTimes<outcomeTimes(tt)+post)),1);
            
            % changing raster color based on trial type
            if isequal(outcomeType(tt),1)
                rasCol = rgb('forestgreen');
            elseif isequal(outcomeType(tt),2)
                rasCol = rgb('orangered');
            end
            
            % plotting rasters for conflict (in the least efficient way possible)
            for sp = 1:size(spikes.channel(ch).unit(un).trial(tt).times,1)
                line([spikes.channel(ch).unit(un).trial(tt).times(sp)-pre spikes.channel(ch).unit(un).trial(tt).times(sp)-pre], [sortedTrialIdcs(tt)-(10/20) sortedTrialIdcs(tt)+(10/20)],'linewidth',1, 'color', rasCol)
            end
		end

       % calculating psths
       kernelWidth = 25  ./1000;
       [Rbank,t,Ebank] = psth(spikes.channel(ch).unit(un).trial(outcomeType==1), kernelWidth, 'n', [0 pre+post]);
       [Rpop,t,Epop] = psth(spikes.channel(ch).unit(un).trial(outcomeType==2), kernelWidth, 'n', [0 pre+post]);

	   % timing
       tsec = t-repmat(pre,1,length(t));

		for t2 = 1:nTrials	   
			% single trial firing rates. 
			if isempty(spikes.channel(ch).unit(un).trial(t2).times)
				Rst(t2,:) = zeros(1,length(t));
				% fprintf('\n     zero spikes recorded in trial %d! \n', t2)
			else
				[Rst(t2,:)] = psth(spikes.channel(ch).unit(un).trial(t2), kernelWidth, 'n', [0 pre+post]);
			end
        end % looping over trials

        axis square
        if ch<=8
			title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{1}))
		elseif ch>8 & ch<=16
		    title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{2}))
		elseif ch>16 & ch<=24
		    title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,microLabels{3}))
		elseif ch>24 & ch<=32
		    title(sprintf('chan: %d, unit: %d; %s',inclChans(ch),un,'4th BF???'))
		else
		    error('trying to loop over more channels than recorded electrodes. something went wrong...')
		end
        xlim([-1 2])
        ylabel('trials')
		
        % plotting the firing rates and SEMs.
        subplot(3,2,3)
        hold on
   
   		try 
        % plotting easy psth
        patch([tsec fliplr(tsec)],[Rbank+Ebank fliplr(Rbank-Ebank)], rgb('forestgreen'),'edgecolor','none','facealpha',0.5)
        plot(tsec,Rbank,'color',rgb('forestgreen'),'linewidth',2)
        
        patch([tsec fliplr(tsec)],[Rpop+Epop fliplr(Rpop-Epop)], rgb('orangered'),'edgecolor','none','facealpha',0.5)
        plot(tsec,Rpop,'color',rgb('orangered'),'linewidth',2)
		end

        hold off
        axis tight square
        xlim([-1 2])
        xlabel('time relative to outcome (s)')
        ylabel('firing rate (spks/s)')
%         halfMaximize(un,'left')
        
		%% ~~~ plotting the prospect curves ~~~
		% Summary statistics for each channel in time.
		meanBank = squeeze(mean(Rst(outcomeType==1,tsec>0.25 & tsec<1.25),2));
		errBank = squeeze(std(Rst(outcomeType==1,tsec>0.25 & tsec<1.25),0,2))./sqrt(nBanks);
		meanPop = squeeze(mean(Rst(outcomeType==2,tsec>0.25 & tsec<1.25),2));
		errPop = squeeze(std(Rst(outcomeType==2,tsec>0.25 & tsec<1.25),0,2))./sqrt(nPops);

		% color map for trial identity
		% outcomeCmap = repmat(rgb('forestgreen'),nTrials,1);
		% outcomeCmap(outcomeType==2,:) = repmat(rgb('orangered'),sum(outcomeType==2),1);
		sortedOutcomeCmap = [repmat(rgb('forestgreen'),sum(outcomeType==1),1); repmat(rgb('orangered'),sum(outcomeType==2),1)];

		%% high gamma vs outcome
		subplot(3,2,2)
		% first fit the data
		lineBuffer = 20;
		lmPop = fitlm(-pointsPerTrial(outcomeType==2),meanPop,'y ~ 1 + x1');
		qmPop = fitlm(-pointsPerTrial(outcomeType==2),meanPop,'y ~ 1 + x1^2');
		fPopLine = polyval(fliplr(lmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);
		fPopQuad = polyval(fliplr(qmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);

		lmBank = fitlm(pointsPerTrial(outcomeType==1),meanBank,'y ~ 1 + x1');
   		qmBank = fitlm(pointsPerTrial(outcomeType==1),meanBank,'y ~ 1 + x1^2');
		fBankLine = polyval(fliplr(lmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);
		fBankQuad = polyval(fliplr(qmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);

		% then plot scatters and fitted lines...
		hold on
		scatter([pointsPerTrial(outcomeType==1) -pointsPerTrial(outcomeType==2)],[meanBank' meanPop'],5,sortedOutcomeCmap,'filled')
		% plot linear fits
		plot(min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer,fPopLine,'color',rgb('orangered'))
		plot(min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer,fBankLine,'color',rgb('forestgreen'))

		% plot quadratic fits
		plot(min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer,fPopQuad,'color',rgb('orangered'),'linestyle','--')
		plot(min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer,fBankQuad,'color',rgb('forestgreen'),'linestyle','--')

		hold off

		% plot three deets
		title('per trial high gamma vs. points gained-lost')
		axis square tight

		%% statistics
		subplot(3,2,4)
		hold on
		betterBoxplot(1,lmPop.Residuals.Raw.^2,rgb('orangered'),4)
		betterBoxplot(2,lmBank.Residuals.Raw.^2,rgb('forestgreen'),4)
		hold off

		title('Squared residuals from each linear model')
		xlim([0.5 2.5])
		axis square

		% printing pop stats
		subplot(3,2,5)
		hold on
		text(0,1,'Linear Pop Model (solid):','fontweight','bold')
		text(0,0.8,sprintf('intercept: %.2f',lmPop.Coefficients.Estimate(1)))
		text(0,0.6,sprintf('slope: %.2f',lmPop.Coefficients.Estimate(2)))
		if lmPop.Coefficients.pValue(2)<0.05
		text(0,0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',lmPop.DFE,lmPop.Coefficients.tStat(2),lmPop.Coefficients.pValue(2)))
		else
		text(0,0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',lmPop.DFE,lmPop.Coefficients.tStat(2),lmPop.Coefficients.pValue(2)))
		end
		text(0,0.2,'Quadratic Pop Model (dashed):','fontweight','bold')
		text(0,0,sprintf('intercept: %.2f',qmPop.Coefficients.Estimate(1)))
		text(0,-0.2,sprintf('slope: %.2f',qmPop.Coefficients.Estimate(2)))
		if qmPop.Coefficients.pValue(2)<0.05
		text(0,-0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',qmPop.DFE,qmPop.Coefficients.tStat(2),qmPop.Coefficients.pValue(2)))
		else
		text(0,-0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',qmPop.DFE,qmPop.Coefficients.tStat(2),qmPop.Coefficients.pValue(2)))
		end
		axis off
		hold off

		% printing bank stats
		subplot(3,2,6)
		hold on
		text(0,1,' Linear Bank Model (solid):','fontweight','bold')
		text(0,0.8,sprintf('intercept: %.2f',lmBank.Coefficients.Estimate(1)))
		text(0,0.6,sprintf('slope: %.2f',lmBank.Coefficients.Estimate(2)))
		if lmBank.Coefficients.pValue(2)<0.05
		text(0,0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',lmBank.DFE,lmBank.Coefficients.tStat(2),lmBank.Coefficients.pValue(2)))
		else
		text(0,0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',lmBank.DFE,lmBank.Coefficients.tStat(2),lmBank.Coefficients.pValue(2)))
		end
		text(0,0.2,'Quadratic Bank Model (dashed):','fontweight','bold')
		text(0,0,sprintf('intercept: %.2f',qmBank.Coefficients.Estimate(1)))
		text(0,-0.2,sprintf('slope: %.2f',qmBank.Coefficients.Estimate(2)))
		if qmBank.Coefficients.pValue(2)<0.05
		text(0,-0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',qmBank.DFE,qmBank.Coefficients.tStat(2),qmBank.Coefficients.pValue(2)))
		else
		text(0,-0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',qmBank.DFE,qmBank.Coefficients.tStat(2),qmBank.Coefficients.pValue(2)))
		end
		axis off
		hold off

		% colormap(turbo)
		halfMaximize(un,'left')

        % saving figures
        saveDir = sprintf('~/Dropbox/BART_firing_prospect/');
        if ~isempty(Rbank)
            %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
            saveas(un,fullfile(saveDir,sprintf('pt%s_ch%d_un%d_banksvspops_firingRate_prospect.pdf',ptID,ch,un)))
        end
        close(un)
        
    end % looping over units
end % looping over channels



















%
%
%
% %% loading kilosorted spikes
% spikesDir = '/media/user1/data4TB/data/spikes_BART_202002';
% spikeStruct = loadKSdir(spikesDir);
% clusterIDs = unique(spikeStruct.clu);
% nUnits = length(clusterIDs);
%
%
% %% spike parameters
% binSize = 10./1000;
% kernelWidth = 50./1000;
%
%
% alignSpots = outcomeTimes;
% % epoch Firing Rates
% for un = 1:nUnits
%     thisClust = spikeStruct.clu;
%     unitTimes = (spikeStruct.st(spikeStruct.clu==clusterIDs(un))');
%
%     % compute everything
%     [psth, bins, rasterX, rasterY, ~, ba] = psthAndBA(unitTimes, alignSpots, [-pre post], binSize);
%
%     % PSTH smoothing filter
%     gw = gausswin(round(kernelWidth*600),3); % kernelwidth[*4]?? in units of ms here, because gausswin fuction likes ints.
%     smWin = gw./sum(gw);
%
%     % smooth ba
%     baSm = conv2(smWin,1,ba', 'same')'./binSize;
%
%     % stacking smoothed psths into a tensor
%     softNorm = false;
%     if softNorm % [20181114] this is what churchland does for all of his papes.
%         trialTensor(un,:,:) = baSm./(range(range(baSm))+5);
%     else
%         trialTensor(un,:,:) = baSm;
%     end
%     % [20181114] Better to observe raw firing rates to start.
%     baTensor(un,:,:) = ba;
%
%     % timing vector
%     tSec = linspace(-pre,post,size(trialTensor,3));
%
%     % visualize the data here for each unit
%     figure(un)
%     subplot(2,1,1)
%     plotSpikeRaster(logical(ba),'PlotType','vertline');
%     xlim([200 400])
%     axis square off
%     title(sprintf('unitID: %d',clusterIDs(un)))
%
%     subplot(2,1,2)
%     hold on
%     patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==2,:),2)-(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))'...
%         fliplr(squeeze(mean(trialTensor(un,outcomeType==2,:),2)+(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))')]...
%         ,rgb('orangered'),'facealpha',0.5,'edgecolor','none')
%     plot(tSec,squeeze(mean(trialTensor(un,outcomeType==2,:),2))','color',rgb('orangered'))
%     text(0,2.5,sprintf('N = %d',sum(outcomeType==2)),'color',rgb('orangered'))
%
%     patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==1,:),2)-(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))'...
%         fliplr(squeeze(mean(trialTensor(un,outcomeType==1,:),2)+(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))')]...
%         ,rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
%     plot(tSec,squeeze(mean(trialTensor(un,outcomeType==1,:),2))','color',rgb('forestgreen'))
%     text(0,3.5,sprintf('N = %d',sum(outcomeType==1)),'color',rgb('forestgreen'))
%     %     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
%     hold off
%     axis tight square
%     xlim([-1 2])
%     xlabel('time relative to outcome (s)')
%     ylabel('firing rate (spks/s)')
%     title(sprintf(''))
%     halfMaximize(un,'left')
%
%     % saving figures
%     saveDir = sprintf('~/Figs/BART/%s',ptID);
%     if exist(saveDir,'dir')
%         %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
%         saveas(un,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_firing rates.pdf',ptID,deblank(trodeLabels{ch2}))))
%     end
%
%
% end
%
%
%
