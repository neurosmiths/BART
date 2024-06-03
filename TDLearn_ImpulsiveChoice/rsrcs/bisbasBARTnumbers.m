function [ptArray,allBHV,hazEEG] = bisbasBARTnumbers(ROI)
% BARTNUMBERS outputs overall summary stats for BART
%
% [ptArray,allBHV,hazEEG] = BARTnumbers() will output a arrays of patient names, behavioral stats, and a logical array of which patients have scalp EEG contacts
% 
% including a string as an input argument will search for matches in an array of contact labels and output the number of contacts with the corresponding anatomical label.

% Edited RC20220413 added new patients. 

if nargin<1
	ROI = '';
end

% numbers of patients
ptArray = {'202114','202117','202118','202201', '202207', '202212',...
    '202214', '202215', '202216', '202217', '202302'}';

% 202212b??

% pt removed from analysis due to low numbers of trials: 202906
% r - indicates revision: same patient; different 'trodes.
% u - indicates file with microelectrodes (stopped doing this pretty quickly)

fprintf('\nnumber of subjects: %d\n',length(ptArray));

if gt(nargout,2)
    for p = 1:length(ptArray)
        allBHV(p) = BARTbehavior(ptArray{p},false);
        
        
        [labels{p},isECoG{p},isEEG{p},isECG{p},anatomicalLocs,~] = ptTrodesBART(ptArray{p});

 		if p==1
 			allBrainAreas = anatomicalLocs;
 		else
 			if size(anatomicalLocs,2)>1
 				allBrainAreas = cat(1,allBrainAreas,anatomicalLocs{:,1});
 			else
 				allBrainAreas = cat(1,allBrainAreas,anatomicalLocs);
 			end
 		end
 
        hazEEG(p) = logical(sum(isEEG{p}));
        
        fprintf('\nimpulsivity metrics for subject %s : %.2f\n',deblank(ptArray{p}),allBHV(p).impulsivityKLD);
    end

fprintf('\nmean +/- std number of trials for %d patients: %.2f +/- %.2f \n',length(ptArray),mean([allBHV.totalTrials]), std([allBHV.totalTrials]));

fprintf('\nmean +/- std accuracy for %d patients: %.2f +/- %.2f \n',length(ptArray),mean([allBHV.accuracyTot]), std([allBHV.accuracyTot]));

fprintf('\ntotal number of electrodes: %d\n',sum(cell2mat(isECoG')))

fprintf('\nmean +/- std number of electrodes: %.2f +/- %.2f \n',sum(cell2mat(isECoG'))./length(ptArray),std(cellfun(@sum,isECoG)))

BAOI = allBrainAreas(cellfun('isclass',allBrainAreas,'char'));

fprintf('\n%d total contacts matching "%s"',sum(contains(BAOI,ROI)),ROI)

fprintf('\n%d out of %d patients with EEG.',sum(hazEEG),length(ptArray))

end


