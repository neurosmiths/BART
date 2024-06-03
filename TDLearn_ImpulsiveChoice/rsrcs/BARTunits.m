function [BARTstats] = BARTunits(ptID,nevFile)
% BARTUNITS analyzes and visualizes firing rates for the BART task.
%
%   [BARTstats] = BARTunits(ptID,nevFile) analyzes firing rate data for the
%   patient specified in the string ptID using the data in nevFile.
%
%   Currently only supports tab delimited text files exported from offline
%   sorter.
%

% nevFile = '/home/elliot/Code/MATLAB/BART/BNST/BNST_recording_OR/20170713-20170713-114259-114259.nev'
% ptID = 'BNST1'

% author: EHS20170713

% just throwing in these files here, now. 
ptID = 'BNST';
nevFile = '/home/user1/code/matlab/BART/BNST/BNST_recording_OR/20170713-20170713-114259-114259.mat'

% timing details.
pre = 3;
post = 3;

% load and define triggers from nevFle
%NEV = openNEV(nevFile)
load(nevFile) % in case nevFile is a mat file
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% loading spike-sorted data and aligning on NEV triggers
nChans = 1;
if ~isequal(nChans,0)
    for ch = 2
        disp(sprintf('organizing spike times for channel %d.',ch))
        % just support for tab delimited text files from offline sorter
        sortedFile = sprintf('%s-SortedEHS-Channel%d.txt',nevFile(1:end-4),ch);
        dtz = importdata(sortedFile);
        UnitTimestamp = dtz(:,2:3);
        
        % banks and pops
        alignName = 'banksandpops';
        bankTimes = trigTimes(trigs==25);
        popTimes = trigTimes(trigs==26);
        nBanks = length(bankTimes);
        nPops = length(popTimes);
        outcomeType = [ones(1,nBanks) 2*ones(1,nPops)]; % 1 = bank, 2 = pop
        nTrials = length(outcomeType);
        
        % looping over units.
        nUnits = length(unique(UnitTimestamp(:,1)));
        for un = 1:nUnits
            unitTimes = UnitTimestamp(UnitTimestamp(:,1)==un,2); % in seconds
            % loooping over banks
            for bnk = 1:nBanks
                % putting the data in a structure
                data.channel(ch).unit(un).banks(bnk).times = unitTimes(unitTimes>bankTimes(bnk)-pre & unitTimes<bankTimes(bnk)+post) - repmat(bankTimes(bnk)-pre,length(unitTimes(unitTimes>bankTimes(bnk)-pre & unitTimes<bankTimes(bnk)+post)),1);
            end
            % loooping over Pops
            for pp = 1:nPops
                % putting the data in a structure
                data.channel(ch).unit(un).pops(pp).times = unitTimes(unitTimes>popTimes(pp)-pre & unitTimes<popTimes(pp)+post) - repmat(popTimes(pp)-pre,length(unitTimes(unitTimes>popTimes(pp)-pre & unitTimes<popTimes(pp)+post)),1);
            end
            
            % now plotting per-bank/pop rasters and rates for each unit,.
            fig1 = figure(un)
	    fig1.Renderer = 'Painters';

	    % plotting rasterds....
            ah_ras = plotmultipleaxes(1,1,2,0.08,un);
            hold on
            for tt = 1:nTrials
                switch outcomeType(tt)
                    case 1
                        rasCol = [0 0.5 0];
                        outcomeLabel = 'banked';
                    case 2
                        rasCol = [1 0.5 0];
                        outcomeLabel = 'popped';
                end
                
                % plotting rasters for conflict (in the least efficient way possible)
                if tt<=nBanks
                    for sp = 1:size(data.channel(ch).unit(un).banks(tt).times,1)
                        try
                            line([data.channel(ch).unit(un).banks(tt).times(sp)-pre data.channel(ch).unit(un).banks(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', rasCol)
                        catch
                            line([data.channel(ch).unit(un).banks(tt).times(sp)-pre data.channel(ch).unit(un).banks(tt).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', 'k')
                        end
                    end
                else
                    for sp = 1:size(data.channel(ch).unit(un).pops(tt-nBanks).times,1)
                        try
                            line([data.channel(ch).unit(un).pops(tt-nBanks).times(sp)-pre data.channel(ch).unit(un).pops(tt-nBanks).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', rasCol)
                        catch
                            line([data.channel(ch).unit(un).pops(tt-nBanks).times(sp)-pre data.channel(ch).unit(un).pops(tt-nBanks).times(sp)-pre], [tt-(9/20) tt+(9/20)],'linewidth',2, 'color', 'k')
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
                ylabel('Trials','fontsize', 16)
                set(gca, 'linewidth', 2, 'fontsize', 16);
                
            end
            hold off
            
            % calculating psths
            kernelWidth = 50  ./1000;
            [Rbank,~,Ebank] = psth(data.channel(ch).unit(un).banks, kernelWidth, 'n', [0 pre+post]);
            [Rpop,t,Epop] = psth(data.channel(ch).unit(un).pops, kernelWidth, 'n', [0 pre+post]);
            tsec = t-repmat(pre,1,length(t));
            
            figure(un)
            ah_ras = plotmultipleaxes(2,1,2,0.08,un);
            hold on
            patch([tsec fliplr(tsec)],[Rbank+Ebank fliplr(Rbank-Ebank)], [0 0.5 0],'edgecolor','none','facealpha',0.5)
            plot(tsec,Rbank,'color',[0 0.5 0],'linewidth',2)
            patch([tsec fliplr(tsec)],[Rpop+Epop fliplr(Rpop-Epop)], [1 0.5 0],'edgecolor','none','facealpha',0.5)
            plot(tsec,Rpop,'color',[1 0.5 0],'linewidth',2)
            
            % stimulus timing lines
%             line([0 0], [0 nTrials],'linestyle', '--', 'color', 'k')
            % PSTH plot details
            xlim([-(pre-1) (post-1)])
            hold off
            xlabel('Time (seconds)', 'fontsize', 16);
            ylabel('Firing Rate (spikes/second)', 'fontsize', 16);
            set(gca, 'linewidth', 2, 'fontsize', 16);
            
            maximize(un)
            fName = sprintf('%s_channel%d_unit%d_firingRates_%s.pdf',ptID,ch,un,'bankpop')
            saveas(un,fName)
            close(un)
        end
    end
else; disp('no units sorted, eh? Elliot will address this soon...'); end
