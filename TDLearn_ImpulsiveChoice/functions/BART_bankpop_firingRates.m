function [nUnits,nTrials] = BART_countUnits(ptID)
% BART_COUNTUNITS provides unit and trial counts for a particular subject
%
%   [nUnits,nTrials] = BART_countUnits(ptID) provides unit and trial counts from a BART session
%	for the subject [ptID].
% 


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
% [trodeLabels,isECoG,isEEG,isECG,anatomicalLocs,adjacentChanMat] = ptTrodesBART(ptID);

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
TimeRes = NEV.MetaTags.TimeRes;

% task parameters in chronological order..
respTimes = trigTimes(trigs==24);
outcomeTimes = trigTimes(trigs==25 | trigs==26);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24; % 1 = bank, 2 = pop
[~,sortedTrialIdcs] = sort(outcomeType);
nTrials = length(outcomeType);

% standard 3-D [chan, unit, timestamp (seconds)] matrix. 
ChanUnitTimestamp = [double(NEV.Data.Spikes.Electrode)' double(NEV.Data.Spikes.Unit)' (double(NEV.Data.Spikes.TimeStamp)./TimeRes)'];

% channel deets. 
inclChans = unique(ChanUnitTimestamp(:,1));
microLabels = microLabelsBART(ptID)
inclChans(inclChans-96>length(microLabels)*8) = []; % magic numbers for recording on bank D and number of BF micros.
nChans = length(inclChans);


% 









































