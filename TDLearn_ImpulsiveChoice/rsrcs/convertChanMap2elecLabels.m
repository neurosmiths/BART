function [labels] = convertChanMap2elecLabels(ChanMap,ptID)
% CONVERTCHANMAP2ELECLABELS converts chanMap data to a csv of electrode
%   labels. 
%
%   chanMap is derived from Tyler Davis' Locate Electrodes GUI. 
%
%   ptID is a string that identifies the patient. 
% 



for ch = 1:length(ChanMap.ElecMapRaw)
    try
        elecLabels{ch,1} = ChanMap.LabelMap{ChanMap.ChannelMap1==ch};
    catch
        elecLabels{ch,1} = 'NaC';
    end
end

% converting to character cell.
labels = cellstr(elecLabels);

% converting to table
T = cell2table(elecLabels);

writetable(T,['D:\Data\preProcessed\BART_preprocessed\electrodeLabels\' ptID 'labels.csv'])