function [microLabels,microPts] = microLabelsBART(ptID)

% MICROLABELSBART outputs channel labels for microwires

% author: EHS20200624

microPts = {'202001','202002','202006','202007','202009','202011','202014',...
    '202015','202016','202105','202107'}

if any(contains(microPts,ptID))
	switch ptID
		case {'202001'}
			microLabels = {'left Amygdala','left Anterior Cingulate'};
		case {'202002'}
			microLabels = {'left medial Orbital Gyrus','left Anterior Cingulate','left Amygdala'};
		case {'202006'}
			microLabels = {'left hippocampus','left dorsal Anterior Cingulate'};
		case {'202007'}
			microLabels = {'left subcallosal area (vmPFC)','left anterior hippocampus'};
		case {'202009'}
			microLabels = {'right Gyrus Rectus','right dorsal Anterior Cingulate'};
        case {'202011'}
			microLabels = {'right Gyrus Rectus','right parahippocampal Gyrus'};
        case {'202014'}
			microLabels = {'right OFC','right hippocampus'};
		case {'202015'}
			microLabels = {'right OFC','right hippocampus'};	
		case {'202016'}
			microLabels = {'right OFC','right hippocampus'};
        case {'202105'} 
            microLabels = {'left OFC','right hippocampus'};
		case {'202107'} 
			microLabels = {'left Subgenual Cingulate','left Anterior Cingulate'};
        case {'202110'}
            microLabels = {'left OFC','left Subgenual Cingulate'};
        case {'202314'}
            microLabels = {'left OFC','left Anterior Cingualte', 'right Amygdala'};
	end
else
	fprintf('\nThis patient may not have had micros...\n')
	microLabels = {};
end
