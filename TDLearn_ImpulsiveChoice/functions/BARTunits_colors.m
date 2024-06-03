function [BARTstats] = BARTunits_colors%(ptID,nevFile)
% BARTUNITS analyzes and visualizes firing rates for the BART task.
%
%   [BARTstats] = BARTunits(ptID,nevFile) analyzes firing rate data for the
%   patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.
%

% author: EHS20170713


% input args
nevFile = '/media/user1/data4TB/data/BART/testData/datafile003.nev'
ptID = 'Elliot'

% timing details.
pre = 3;
post = 3;

% load and define triggers from nevFle
NEV = openNEV(nevFile);
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% balloon colors
alignName = 'colors';
[outcomeType,sortIDx] = sort(trigs(trigs<20));
trialStarts = trigTimes(trigs<20);
trialStarts = trialStarts(sortIDx);
nTrials = length(outcomeType);

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% loading spike-sorted data and aligning on NEV triggers
nChans = input('how many unit channels were recorded and sorted?');
if ~isequal(nChans,0)
    for ch = 1:nChans
        disp(sprintf('organizing spike times for channel %d.',ch))
        % just support for tab delimited text files from offline sorter
        sortedFile = sprintf('%s-SortedEHS-Channel%d.txt',nevFile(1:end-4),ch);
        dtz = importdata(sortedFile);
        UnitTimestamp = dtz(:,2:3);
        
        % looping over units.
        nUnits = length(unique(UnitTimestamp(:,1)));
        for un = 1:nUnits
            unitTimes = UnitTimestamp(UnitTimestamp(:,1)==un,2); % in seconds
            % loooping over banks
            for t1 = 1:nTrials
                % putting the data in a structure
                data.channel(ch).unit(un).trials(t1).times = unitTimes(unitTimes>trialStarts(t1)-pre & unitTimes<trialStarts(t1)+post) - repmat(trialStarts(t1)-pre,length(unitTimes(unitTimes>trialStarts(t1)-pre & unitTimes<trialStarts(t1)+post)),1);
            end
            
            % now plotting per-bank/pop rasters and rates for each unit,.
            figure(un)
            ah_ras = plotmultipleaxes(1,1,2,0.08,un);
            hold on
            for tt = 1:nTrials
                for sp = 1:size(data.channel(ch).unit(un).trials(tt).times,1)
                    try
                        line([data.channel(ch).unit(un).trials(tt).times(sp)-pre data.channel(ch).unit(un).trials(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', cMap(outcomeType(tt),:))
                    catch
                        line([data.channel(ch).unit(un).trials(tt).times(sp)-pre data.channel(ch).unit(un).trials(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', cMap(outcomeType(tt)-10,:))
                    end
                end
            end
            
            % stimulus timing lines
            line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
            % raster plot details
            xlim([-(pre-1) (post-1)])
            ylim([0 nTrials])
            str = sprintf('patient %s, Channel %d, Unit %d; aligned on %s',ptID ,ch ,un ,alignName);
            title(str,'fontsize',18);
            ylabel('Trials','fontsize', 14)
            set(gca, 'linewidth', 2, 'fontsize', 12);
            hold off
            
            % calculating psths
            kernelWidth = 50  ./1000;
            [Ryellow,~,Eyellow] = psth(data.channel(ch).unit(un).trials(outcomeType==1), kernelWidth, 'n', [0 pre+post]);
            [Rorange,~,Eorange] = psth(data.channel(ch).unit(un).trials(outcomeType==2), kernelWidth, 'n', [0 pre+post]);
            [Rred,~,Ered] = psth(data.channel(ch).unit(un).trials(outcomeType==3), kernelWidth, 'n', [0 pre+post]);
            [RyellowCTRL,~,EyellowCTRL] = psth(data.channel(ch).unit(un).trials(outcomeType==11), kernelWidth, 'n', [0 pre+post]);
            [RorangeCTRL,~,EorangeCTRL] = psth(data.channel(ch).unit(un).trials(outcomeType==12), kernelWidth, 'n', [0 pre+post]);
            [RredCTRL,~,EredCTRL] = psth(data.channel(ch).unit(un).trials(outcomeType==13), kernelWidth, 'n', [0 pre+post]);
            [RgrayCTRL,t,EgrayCTRL] = psth(data.channel(ch).unit(un).trials(outcomeType==14), kernelWidth, 'n', [0 pre+post]);
            
            tsec = t-repmat(pre,1,length(t));
            
            figure(un)
            ah_ras = plotmultipleaxes(2,1,2,0.08,un);
            hold on
%             patch([tsec fliplr(tsec)],[Ryellow+Eyellow fliplr(Ryellow-Eyellow)], cMap(1,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[Rorange+Eorange fliplr(Rorange-Eorange)], cMap(2,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[Rred+Ered fliplr(Rred-Ered)], cMap(3,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[RyellowCTRL+EyellowCTRL fliplr(RyellowCTRL-EyellowCTRL)], cMap(11-10,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[RorangeCTRL+EorangeCTRL fliplr(RorangeCTRL-EorangeCTRL)], cMap(12-10,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[RredCTRL+EredCTRL fliplr(RredCTRL-EredCTRL)], cMap(13-10,:),'edgecolor','none','facealpha',0.5)
%             patch([tsec fliplr(tsec)],[RgrayCTRL+EgrayCTRL fliplr(RgrayCTRL-EgrayCTRL)], cMap(14-10,:),'edgecolor','none','facealpha',0.5)
            
            plot(tsec,Ryellow,'color', cMap(1,:),'linewidth',2)
            plot(tsec,Rorange,'color', cMap(2,:),'linewidth',2)
            plot(tsec,Rred,'color', cMap(3,:),'linewidth',2)
            plot(tsec,RyellowCTRL,'color', cMap(11-10,:),'linewidth',2,'linestyle','--')
            plot(tsec,RorangeCTRL,'color', cMap(12-10,:),'linewidth',2,'linestyle','--')
            plot(tsec,RredCTRL,'color', cMap(13-10,:),'linewidth',2,'linestyle','--')
            plot(tsec,RgrayCTRL,'color', cMap(14-10,:),'linewidth',2,'linestyle','--')

            % stimulus timing lines
            %             line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
            % PSTH plot details
            xlim([-(pre-1) (post-1)])
            hold off
            xlabel('time relative to balloon presentation (s)', 'fontsize', 14);
            ylabel('firing rate (spikes/s)', 'fontsize', 14);
            set(gca, 'linewidth', 2, 'fontsize', 12);
            
            fDir = '~/Dropbox/Figs/BNST_BART/';
            fName = sprintf('%s_channel%d_unit%d_firingRates_balloonColor_%s.pdf',ptID,ch,un,alignName);
            print(un,[fDir fName],'-dpdf','-fillpage')
            close(un)
        end
    end
else; disp('no units sorted, eh? Elliot will address this soon...'); end