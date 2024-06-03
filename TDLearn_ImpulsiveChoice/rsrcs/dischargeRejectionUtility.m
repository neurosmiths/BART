function [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionUtility(ptID,LFPmat,plotFlag)
%DISCHARGREJECTIONUTILITY reject trials with epileptiform discharges.
%
%   [rejectTrialsperChannel,trialsAcrossChannels] = dischargeRejectionTool
%       (ptID,LFPmat) finds outliers based on their maximum and mimimum values.
%
%       LFPmat should be organized as [channels X samples X trials]
%
%   optional input argument: plotFlag == 0 will not plot rejected trials.
%

% author: ElliotHSmith (https://github.com/elliothsmith/MSIT-analysis)

if plotFlag
    % figure window.
    ah = figure('Color',[0 0 0]);
    hold on
end
try
    [trodeLabels,isECoG,~,~,anatomicalLocs] = ptTrodesBART(ptID);
catch
    fprintf('\n unable to get electrode information for this pt')
end

% data 
nChans = size(LFPmat,1);
time = 1:size(LFPmat,2);
nTrials = size(LFPmat,3);

for ch = 1:nChans
    fprintf('\nremoving interictal discharges for patient %s, channel %s',ptID,deblank(trodeLabels{ch}))
    % rejecting trials using the standard outlier function
    [rejectTrialsperChannel(ch).rejectTheseTrials] = outliers(range(squeeze(LFPmat(ch,:,:))));
    
    % saving all of the bad trials
    if isequal(ch,1)
        trialsAcrossChannels = find(outliers(range(squeeze(LFPmat(ch,:,:)))));
    else
        trialsAcrossChannels = cat(2,trialsAcrossChannels,find(outliers(range(squeeze(LFPmat(ch,:,:))))));
    end
    
    % plotting each trial, like a dunce...
    for tt = 1:nTrials
        if plotFlag
            %% [20160622] plotting LFP from each trial brushing to reject
            if rejectTrialsperChannel(ch).rejectTheseTrials
                plot(time,LFPmat(ch,:,tt)+((tt-1)*1000),'color',rgb('orangered'));
            else
                plot(time,LFPmat(ch,:,tt)+((tt-1)*1000),'color',rgb('white'));
            end
            axis tight off
        end
    end
%     saveas(gcf,sprintf('interictal discharges for patient %s, channel %s',ptID,deblank(trodeLabels{ch}))
end % looping over lfp channels

trialsAcrossChannels = unique(trialsAcrossChannels);
