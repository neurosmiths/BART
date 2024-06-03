function [labels] = convertChanMap2elecLabels_newGUI(fname,ptID)
% CONVERTCHANMAP2ELECLABELS converts chanMap data to a csv of electrode
%   labels. 
%
%   fname refers to the appropriate ChannelMap.m is derived from Tyler 
%       Davis' Locate Electrodes GUI. 
%
%   ptID is a string that identifies the patient. 
% 

% loading file
load(fname)

% looping over and assigning electrodes. 
for ch = 1:length(ElecXYZProj)
    if str2double(ptID)>=202015
        try
            elecLabels{ch,1} = LabelMap{ChannelMap1==ch};
        catch
            elecLabels{ch,1} = 'NaC';
        end
    else
        try
            elecLabels{ch,1} = LabelMap{ChannelMap1==ch};
        catch
            elecLabels{ch,1} = 'NaC';
        end
    end
end

% converting to character cell.
labels = cellstr(elecLabels);

% converting to table
T = cell2table(elecLabels);

writetable(T,fullfile('D:\','Data','preProcessed','BART_preprocessed','electrodeLabels_wEEG',sprintf('%slabels.csv',ptID)))
