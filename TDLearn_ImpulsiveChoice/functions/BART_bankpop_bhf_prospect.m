function [] = BART_bankpop_bhf_prospect(ptID)   % BARTstats
% BART_BANKPOP_BHF analyzes and visualizes LFP for the BART task.
%
%   [BARTstats] = BART_bankpop_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.

% author: EHS20181005


%  ptID = '201902r';


nevList = dir(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/*.nev',ptID));
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
nChans = length(trodeLabels(isECoG));
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% timing parameters.
pre = 3;
post = 3;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
pointsPerTrial = [data.points];
if length(balloonTimes)>length(outcomeTimes)
    balloonTimes(end) = [];
end
if length(inflateTimes)>length(outcomeTimes)
    inflateTimes(end) = [];
end
if length(pointsPerTrial)>length(outcomeTimes)
    pointsPerTrial(end) = [];
end
nTrials = length(outcomeType);

% HG filter
[b,a] = butter(4,[70 200]/(Fs/2));

% epoching data
LFPmat = zeros(nChans,Fs*(pre+post)+1,nTrials);
for tt = 1:nTrials
    % epoch the data here [channels X samples X trials]
    LFPmat(:,:,tt) = NSX.Data(1:nChans,floor(Fs*outcomeTimes(tt))-Fs*pre:floor(Fs*outcomeTimes(tt))+Fs*post);
    % do spectral claculations here
    for ch = 1:nChans
        HGmat(ch,:,tt) = abs(hilbert(filtfilt(b,a,double(LFPmat(ch,:,tt)))));
    end
end


%% TODO:: remove HG outliers


% how much smoothing in time? (50 - 100  ms is usually good for BHF) 
smoothType = 'movmean';
smoothFactor = Fs./5;
% Summary statistics for each channel in time. 
nBanks = sum(outcomeType==1); 
nPops = sum(outcomeType==2);
meanBank = smoothdata(squeeze(mean(HGmat(:,:,outcomeType==1),3)),2,smoothType,smoothFactor)';
errBank = smoothdata(squeeze(std(HGmat(:,:,outcomeType==1),0,3)),2,smoothType,smoothFactor)'./sqrt(nBanks);
meanPop = smoothdata(squeeze(mean(HGmat(:,:,outcomeType==2),3)),2,smoothType,smoothFactor)';
errPop = smoothdata(squeeze(std(HGmat(:,:,outcomeType==2),0,3)),2,smoothType,smoothFactor)'./sqrt(nPops);

% color map for trial identity
% outcomeCmap = repmat(rgb('forestgreen'),nTrials,1);
% outcomeCmap(outcomeType==2,:) = repmat(rgb('orangered'),sum(outcomeType==2),1);
sortedOutcomeCmap = [repmat(rgb('forestgreen'),sum(outcomeType==1),1); repmat(rgb('orangered'),sum(outcomeType==2),1)];

%% TODO:: calculate reversal point. 


% for visualizing the data.  % plotting bank/pop responses
for ch2 = 1:nChans
    
    % visualize the data here for each channel
    if ishandle(ch2); close(ch2); end
    figure(ch2)
    %% plotting average high gamma responses
    subplot(3,2,1)
    hold on
    patch([tSec fliplr(tSec)],[meanPop(:,ch2)'-errPop(:,ch2)' fliplr(meanPop(:,ch2)'+errPop(:,ch2)')],rgb('orangered'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,meanPop(:,ch2)','color',rgb('orangered'))
    text(-2,max(meanPop(:,ch2)'+errPop(:,ch2)')-1,sprintf('N = %d',uint16(nPops)),'color',rgb('orangered'))
    patch([tSec fliplr(tSec)],[meanBank(:,ch2)'-errBank(:,ch2)' fliplr(meanBank(:,ch2)'+errBank(:,ch2)')],rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
    plot(tSec,meanBank(:,ch2)','color',rgb('forestgreen'))
    text(-2,max(meanPop(:,ch2)'+errPop(:,ch2)'),sprintf('N = %d',uint16(nBanks)),'color',rgb('forestgreen'))
%     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
    hold off
    
    % plot one deets
    axis tight square
    xlim([-2 2])
    xlabel('time relative to outcome (s)')
    ylabel('broadband high frequency power')
    title([deblank(trodeLabels{ch2}) ' -- ' deblank(anatomicalLocs{ch2})])
    
    
    %% plotting all trials high gamma
    subplot(3,2,2)
    [~,sortedIdcs] = sort(outcomeType);
    hold on
    imagesc(tSec,1:nTrials,smoothdata(squeeze(HGmat(ch2,:,sortedIdcs)),smoothType,smoothFactor)')
    scatter(repmat(2,1,nTrials),1:nTrials,8,sortedOutcomeCmap,'filled','s')
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
    
    lmPop = fitlm(-pointsPerTrial(outcomeType==2),squeeze(mean(HGmat(ch2,tSec>0.25 & tSec<1.25,outcomeType==2),2))','y ~ 1 + x1');
%     qmPop = fitlm(-pointsPerTrial(outcomeType==2),squeeze(mean(HGmat(ch2,tSec>0.25 & tSec<1.25,outcomeType==2),2))','y ~ 1 + x1^2');    
    fPopLine = polyval(fliplr(lmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);
%     fPopQuad = polyval(fliplr(qmPop.Coefficients.Estimate'),min(-pointsPerTrial(outcomeType==2))-lineBuffer:max(-pointsPerTrial(outcomeType==2))+lineBuffer);
    
    lmBank = fitlm(pointsPerTrial(outcomeType==1),squeeze(mean(HGmat(ch2,tSec>0.5 & tSec<2,outcomeType==1),2))','y ~ 1 + x1');
%     qmBank = fitlm(pointsPerTrial(outcomeType==1),squeeze(mean(HGmat(ch2,tSec>0.5 & tSec<2,outcomeType==1),2))','y ~ 1 + x1^2');
    fBankLine = polyval(fliplr(lmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);
%     fBankQuad = polyval(fliplr(qmBank.Coefficients.Estimate'),min(pointsPerTrial(outcomeType==1))-lineBuffer:max(pointsPerTrial(outcomeType==1))+lineBuffer);
    
    
    % then plot scatters and fitted lines...
    hold on
    scatter([pointsPerTrial(outcomeType==1) -pointsPerTrial(outcomeType==2)],...
        [squeeze(mean(HGmat(ch2,tSec>0.5 & tSec<2,outcomeType==1),2))' squeeze(mean(HGmat(ch2,tSec>0.25 & tSec<1.25,outcomeType==2),2))'],...
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

    % printing bank stats
    subplot(3,2,6)
    hold on
    text(0,1,'Bank Model:','fontweight','bold')
    text(0,0.8,sprintf('intercept: %.2f',lmBank.Coefficients.Estimate(1)))
    text(0,0.6,sprintf('slope: %.2f',lmBank.Coefficients.Estimate(2)))
    if lmBank.Coefficients.pValue(2)<0.05
        text(0,0.4,sprintf('SIGNIFICANT model. t(%d) = %.2f, p = %.2f',lmBank.DFE,lmBank.Coefficients.tStat(2),lmBank.Coefficients.pValue(2)))
    else
        text(0,0.4,sprintf('insignficant model. t(%d) = %.2f, p = %.2f',lmBank.DFE,lmBank.Coefficients.tStat(2),lmBank.Coefficients.pValue(2)))
    end
    axis off
    hold off
    
    
    % saving RPE statistics in a structure
    prospectStats(ch2).trodeLabel = deblank(trodeLabels{ch2});
    prospectStats(ch2).PopLinearModel = lmPop;
    prospectStats(ch2).BankLinearModel = lmBank;
        
    
    % saving figures
    saveDir = sprintf('/media/user1/data4TB/Figs/BART/%s/',ptID);
    if exist(fullfile(saveDir,'BHFprospect'),'dir')
        saveas(ch2,fullfile(saveDir,'BHFprospect',sprintf('pt%s_%s_banksvspops_BHF_prospect.pdf',ptID,deblank(trodeLabels{ch2}))))
    else
        mkdir(saveDir,'BHFprospect')
        saveas(ch2,fullfile(saveDir,'BHFprospect',sprintf('pt%s_%s_banksvspops_BHF_prospect.pdf',ptID,deblank(trodeLabels{ch2}))))
    end
    
    close(ch2)

    %% TODO:: do this with theta...
end

% saving stats
save(fullfile(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID),sprintf('pt%s_banksvspops_BHF_prospectStats.mat',ptID)),'prospectStats')


end






%% old code

% for fitting lines using QR decomposition - more efficient, but less
% info/utility wiht regard to understanding model fits. 
%     [pPop,sPop] = polyfit(-pointsPerTrial(outcomeType==2),squeeze(mean(HGmat(ch2,tSec>0.25 & tSec<1.25,outcomeType==2),2))',1);


