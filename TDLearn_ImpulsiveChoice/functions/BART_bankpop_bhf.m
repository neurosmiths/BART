function [statsFile] = BART_bankpop_bhf(ptID)   % BARTstats
% BART_BANKPOP_BHF analyzes and visualizes LFP for the BART task.
%
%   [statsFile] = BART_bankpop_bhf(ptID,nevFile) analyzes LFP data for
%   the patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.

% author: EHS20181005
% Edited: RC20220215

ptID = '202114'

nevList = dir(['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\*.nev']);
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

% load neural data
[nevPath,nevName,nevExt] = fileparts(nevFile);
NSX = openNSx(fullfile(nevPath,[nevName '.ns2']));

% data parameters
nChans = length(trodeLabels(isECoG));
nSamps = size(NSX.Data,2);
Fs = NSX.MetaTags.SamplingFreq;

% timing parameters.
pre = 2;
post = 4;
tSec = linspace(-pre,post,Fs*(pre+post)+1);

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
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

% how much smoothing in time? (50 - 100  ms is usually good for BHF)
smoothFactor = Fs./10;
% Summary statistics for each channel in time.
nBanks = sum(outcomeType==1);
nPops = sum(outcomeType==2);
meanBank = smoothdata(squeeze(mean(HGmat(:,:,outcomeType==1),3)),2,'movmean',smoothFactor)';
errBank = smoothdata(squeeze(std(HGmat(:,:,outcomeType==1),0,3)),2,'movmean',smoothFactor)'./sqrt(nBanks);
meanPop = smoothdata(squeeze(mean(HGmat(:,:,outcomeType==2),3)),2,'movmean',smoothFactor)';
errPop = smoothdata(squeeze(std(HGmat(:,:,outcomeType==2),0,3)),2,'movmean',smoothFactor)'./sqrt(nPops);


permStats = true;
if permStats
    % permutation tests between conditions for each channel
    if ~exist(['\' nevList.folder '\BHFclusterStats_\' ptID '\.mat'],'file')
        tic
        [pval,t_orig,clust_info,~,est_alpha] = clust_perm2(HGmat(:,:,outcomeType==1),HGmat(:,:,outcomeType==2),adjacentChanMat, 500);
        A = toc;
        fprintf('\ncluster-based permutation stats took %d seconds',A)
        save(['\' nevList.folder '\BHFclusterStats_\' ptID '\.mat'],'pval','t_orig','clust_info','est_alpha')

        % for visualizing the data.  % plotting bank/pop responses
        for ch2 = 1:nChans

            % visualize the data here for each channel
            if ishandle(ch2); close(ch2); end
            figure(ch2)
            hold on
            patch([tSec fliplr(tSec)],[meanPop(:,ch2)'-errPop(:,ch2)' fliplr(meanPop(:,ch2)'+errPop(:,ch2)')],rgb('orangered'),'facealpha',0.5,'edgecolor','none')
            plot(tSec,meanPop(:,ch2)','color',rgb('orangered'))
            text(0,2.5,sprintf('N = %d',uint16(nPops)),'color',rgb('orangered'))
            patch([tSec fliplr(tSec)],[meanBank(:,ch2)'-errBank(:,ch2)' fliplr(meanBank(:,ch2)'+errBank(:,ch2)')],rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
            plot(tSec,meanBank(:,ch2)','color',rgb('forestgreen'))
            text(0,3.5,sprintf('N = %d',uint16(nBanks)),'color',rgb('forestgreen'))
            %     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
            hold off
            axis tight square
            xlim([-2 2])
            xlabel('time relative to outcome (s)')
            ylabel('broadband high frequency power')
            title([deblank(trodeLabels{ch2}) ' -- ' deblank(anatomicalLocs{ch2})])

            % saving figures
            saveDir = sprintf('D:\Data\Rhiannon\BART_RLDM_outputs\%s\Figs',ptID);
            if exist(saveDir,'dir')
                %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
                saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
            else
                mkdir(saveDir)
                saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
            end

        end
    else
        fprintf('stats have already been done, and figures have been made. Aren"t you lucky!')
        %load(sprintf('%s/BHFclusterStats_%s.mat',nevList.folder,ptID),'pval','t_orig','clust_info','est_alpha')


    end

    statsFile = sprintf('%s\BHFclusterStats_%s.mat',nevList.folder,ptID)
else
    % for visualizing the data.  % plotting bank/pop responses
    for ch2 = 1:nChans

        % visualize the data here for each channel
        if ishandle(ch2); close(ch2); end
        figure(ch2)
        hold on
        patch([tSec fliplr(tSec)],[meanPop(:,ch2)'-errPop(:,ch2)' fliplr(meanPop(:,ch2)'+errPop(:,ch2)')],rgb('orangered'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanPop(:,ch2)','color',rgb('orangered'))
        text(0,2.5,sprintf('N = %d',uint16(nPops)),'color',rgb('orangered'))
        patch([tSec fliplr(tSec)],[meanBank(:,ch2)'-errBank(:,ch2)' fliplr(meanBank(:,ch2)'+errBank(:,ch2)')],rgb('forestgreen'),'facealpha',0.5,'edgecolor','none')
        plot(tSec,meanBank(:,ch2)','color',rgb('forestgreen'))
        text(0,3.5,sprintf('N = %d',uint16(nBanks)),'color',rgb('forestgreen'))
        %     plot(tSec(logical(pval(ch2,:)<1.9)),pval(logical(pval(ch2,:)<1.9)),'color',rgb('dimgray'),'linewidth',5)
        hold off
        axis tight square
        xlim([-2 2])
        xlabel('time relative to outcome (s)')
        ylabel('broadband high frequency power')
        title([deblank(trodeLabels{ch2}) ' -- ' deblank(anatomicalLocs{ch2})])

        % saving figures
        saveDir = sprintf('D:\Data\Rhiannon\BART_RLDM_outputs\%s\Figs',ptID);
        if exist(saveDir,'dir')
            %         save(sprintf('%s/%s_bankpop_permutationTests.mat',saveDir,ptID),'p')
            saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        else
            mkdir(saveDir)
            saveas(ch2,fullfile(saveDir,sprintf('pt%s_%s_banksvspops_BHF.pdf',ptID,deblank(trodeLabels{ch2}))))
        end

    end

end
