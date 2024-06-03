function [] = BART_bankpop_firingRates_KS(ptID)   % BARTstats
% BART_BANKPOP_FIRINGRATES analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_bankpop_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.

% author: EHS20181005
ptID = '202002';

nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
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

% timing parameters.
pre = 3;
post = 3;

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
nTrials = length(outcomeType);


%% loading kilosorted spikes
spikesDir = '/media/user1/data4TB/data/spikes_BART_202002';
spikeStruct = loadKSdir(spikesDir);
clusterIDs = unique(spikeStruct.clu);
nUnits = length(clusterIDs);


%% spike parameters
binSize = 10./1000;
kernelWidth = 50./1000;


alignSpots = outcomeTimes;
% epoch Firing Rates
for un = 1:nUnits
    thisClust = spikeStruct.clu;
    unitTimes = (spikeStruct.st(spikeStruct.clu==clusterIDs(un))');
    
    % compute everything
    [psth, bins, rasterX, rasterY, ~, ba] = psthAndBA(unitTimes, alignSpots, [-pre post], binSize);
    
    % PSTH smoothing filter
    gw = gausswin(round(kernelWidth*600),3); % kernelwidth[*4]?? in units of ms here, because gausswin fuction likes ints.
    smWin = gw./sum(gw);
    
    % smooth ba
    baSm = conv2(smWin,1,ba', 'same')'./binSize;
    
    % stacking smoothed psths into a tensor
    softNorm = false;
    if softNorm % [20181114] this is what churchland does for all of his papes.
        trialTensor(un,:,:) = baSm./(range(range(baSm))+5);
    else
        trialTensor(un,:,:) = baSm;
    end
    % [20181114] Better to observe raw firing rates to start.
    baTensor(un,:,:) = ba;
    
    % timing vector
    tSec = linspace(-pre,post,size(trialTensor,3));
        
    % visualize the data here for each unit
    figure(un)
    subplot(2,1,1)
    plotSpikeRaster(logical(ba),'PlotType','vertline');
    xlim([200 400])
    axis square off
    title(sprintf('unitID: %d',clusterIDs(un)))
    
    subplot(2,1,2)
    hold on
    patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==2,:),2)-(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))'...
        fliplr(squeeze(mean(trialTensor(un,outcomeType==2,:),2)+(std(trialTensor(un,outcomeType==2,:),[],2)./sqrt(sum(outcomeType==2))))')]...
        ,rgb('orangered'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,squeeze(mean(trialTensor(un,outcomeType==2,:),2))','color',rgb('orangered'))
    text(0,2.5,sprintf('N = %d',sum(outcomeType==2)),'color',rgb('orangered'))
    
    patch([tSec fliplr(tSec)],[squeeze(mean(trialTensor(un,outcomeType==1,:),2)-(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))'...
        fliplr(squeeze(mean(trialTensor(un,outcomeType==1,:),2)+(std(trialTensor(un,outcomeType==1,:),[],2)./sqrt(sum(outcomeType==1))))')]...
        ,rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,squeeze(mean(trialTensor(un,outcomeType==1,:),2))','color',rgb('forestgreen'))
    text(0,3.5,sprintf('N = %d',sum(outcomeType==1)),'color',rgb('forestgreen'))
    %     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
    hold off
    axis tight square
    xlim([-1 2])
    xlabel('time relative to outcome (s)')
    ylabel('firing rate (spks/s)')
    title(sprintf(''))
    halfMaximize(un,'left')
    
    % saving figures
    saveDir = sprintf('~/Figs/BART/%s',ptID);
    if exist(saveDir,'dir')
        %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
        saveas(un,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_firing rates.pdf',ptID,deblank(trodeLabels{ch2}))))
    end
    
    
end



