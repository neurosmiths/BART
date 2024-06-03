% this script visualizes BART anaylses on a flat MNI brain.
close all; clear;

% for faster, non-vectorized graphics on linux...
set(0,'defaultfigurerenderer','painters')


%% ~~~~~~~~~~~~~ FLAT BRAINS ~~~~~~~~~~~~~~~~~

%% organizing [and plotting] surface data.
vflat = []; fFlat = []; fSurf = []; vSurf = []; curv = [];
for hmi = 1:2
    switch hmi
        case 1
            whichHemi = 'left';
            surfPath = '/home/user1/pycortex/filestore/db/MNI/surfaces/pia_lh.gii';
            flatPath = '/home/user1/pycortex/filestore/db/MNI/surfaces/flat_lh.gii';
            curvPath = '/home/user1/code/matlab/brainVisualization/trunk/pycortex_MNI/left.npy';
        case 2
            whichHemi = 'right';
            surfPath = '/home/user1/pycortex/filestore/db/MNI/surfaces/pia_rh.gii';
            flatPath = '/home/user1/pycortex/filestore/db/MNI/surfaces/flat_rh.gii';
            curvPath = '/home/user1/code/matlab/brainVisualization/trunk/pycortex_MNI/right.npy';
    end
    
    % [20190320] load pycortex flat MNI surfaces and freesurfer MNI
    gFlat = gifti(flatPath);
    gSurf = gifti(surfPath);
    
    % load curvature
    curv = readNPY(curvPath);
    
    % translating each hemisphere for dual display
    nVerts = length(gFlat.vertices);
    translationFactor = 180;
    transMat = translationFactor.*[ones(nVerts,1) zeros(nVerts,1)];
    
    % asserting that there are the same number of vertices in the flat and curved surfs
    if ~isequal(nVerts,length(gSurf.vertices))
        error('numbers of vertices do not match between flat and curved MNI surfaces.')
    end
    
    % plotting surfaces
    switch hmi
        case 1
            % left
            vFlat = gFlat.vertices(:,1:2)-transMat;
            vSurf = gSurf.vertices;
            fFlat = gFlat.faces;
            fSurf = gSurf.faces;
            
        case 2
            % right
            vFlat = cat(1,vFlat,gFlat.vertices(:,1:2)+transMat);
            vSurf = cat(1,vSurf,gSurf.vertices);
            fFlat = cat(1,fFlat,gFlat.faces);
            fSurf = cat(1,fSurf,gSurf.faces);
            
    end
end


%% ~~~~~~~~~~~~~ ELECTRODES ~~~~~~~~~~~~~~~~~
% which model and analysis for coloring electrodes.
whichAnalysis = 'FMThetaCorr'; % balloonSizeTrackingGLMresults or bankpopBHF or FMThetaCorr
switch whichAnalysis
    case {'balloonSizeTrackingGLMresults'}
        whichModel = 'each';
    case {'FMThetaCorr'}
        whichModel = '';
    case {'bankpopBHF'}
        whichModel = '';
end
% [20190312] The normalized data will have missing entries for failed
% models.
%% now dealing with sEEG contacts per subject.
[ptArray] = BARTnumbers;
% looping over patients
nPts = length(ptArray);
% Only those patients who have EEG
if strcmp('FMThetaCorr',whichAnalysis)
    pts = [1 2 3 6 8 9 10 11];
else
    pts = 1:nPts;
end
trodeData = struct();
for pt = pts
    % which patient.
    ptID = ptArray{pt};
    
    %% loading [patient-specific] model data and brain surface.
    % loading electrode localization and brain surface data.
    load(sprintf('~/data/BART/BART_EMU/%s/Imaging/Registered/Electrodes.mat',ptID));
    [trodeLabels,isECoG,isEEG] = ptTrodesBART(ptID);
    
    switch whichAnalysis
        case {'balloonSizeTrackingGLMresults'}
            % loading statistics structures from LMMs.
            BARTdir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
            fName = sprintf('%s_%s_%sTrodesModel.mat',ptID,whichAnalysis,whichModel);
            load(fullfile(BARTdir,fName))
            
            % putting data in a struct
            trodeData(pt).patientID = ptID;
            trodeData(pt).trodeLabels = trodeLabels(isECoG);
            trodeData(pt).MNItrodeLocs = ElecXYZMNIProjRaw;
            trodeData(pt).tVals = [trodeLMM.model_t];
            trodeData(pt).pVals = [trodeLMM.model_pCorr];
            trodeData(pt).dfs = [trodeLMM.model_df]';
        case {'FMThetaCorr'}
            % load data from the batch analysis for visualization
            BARTdir = sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data',ptID);
            fName = sprintf('pt%s_Fz_outcome_*_EEGECoGcorrelations_preCuebaselineNorm.mat',ptID);
            trodeFiles = dir(fullfile(BARTdir,fName));
            
            % assert that trodeData is the same length as the numver of
            % electrdoes
            % if strcmp(ptID,'201810')
            % else
            
            if ~isequal(length(trodeLabels(isECoG)),length(trodeFiles))
                error('uh oh! number of trodes in stats does not match that in the electrode lables .csv')
            end
             
            % aggregate data across electrodes
            trodeData(pt).tVals = zeros(1,length(trodeData));
            trodeData(pt).pVals = zeros(1,length(trodeData));
            for el = 1:length(trodeFiles)
                % load data
                load(fullfile(trodeFiles(el).folder,trodeFiles(el).name))
                % find the right trode position in the array
                trodeIdx = contains(trodeLabels(isECoG),EEGECoGtrodeLM.ECoGtrode);
                
                % save data
                trodeData(pt).tVals(trodeIdx) = EEGECoGtrodeLM.trodeLM_BHF.Coefficients{2,3};
                trodeData(pt).pVals(trodeIdx) = EEGECoGtrodeLM.trodeLM_BHF.Coefficients{2,4};
            end
                        
            % putting data in a struct
            trodeData(pt).dfs = EEGECoGtrodeLM.trodeLMtheta.DFE;
            trodeData(pt).patientID = ptID;
            trodeData(pt).trodeLabels = trodeLabels(isECoG);
            trodeData(pt).MNItrodeLocs = ElecXYZMNIProjRaw;
            
    end
    
    % updating user
    fprintf('\nfinshed compiling %s data for patient %d of %d',whichAnalysis,pt,nPts)
    clearvars -except pt nPts ptArray trodeData ...
        whichAnalysis whichModel whichVis ...
        fHandle vFlat vSurf fFlat fSurf
    
end


% data already consolidated for EEG/ECoG comparisons. note: This is based
% on an external .csv
if strcmp('FMThetaCorrPop',whichAnalysis)
    % loading significant contacts.
    A = readtable('/home/user1/Dropbox/EEG_BHFsignificantCorrelations.csv');
    
    % data
    sigData.sigTrodes = A.ECoGtrode;
    sigData.tVals = A.tStat;
    sigData.pVals = A.pVal;
    sigData.rSquared = A.rSquared;
    nSigTrodes = length(A);
end


%% colormap
% making a colormap based o[n stats across all electrodes, all patients
tmpT = horzcat(trodeData.tVals);
tmpT(outliers(tmpT)) = NaN; % removing outliers
posC = sum(tmpT>0);
negC = sum(tmpT<=0);
if posC>negC
    cMapInit = turbo(posC*2);
elseif negC>posC
    cMapInit = turbo(negC*2);
end
cMapInit(cMapInit(:,2)>0.5,2) = 0.4;
tMap = linspace(-max(abs(tmpT)),max(abs(tmpT)),length(cMapInit))';
nTrodes = length(tmpT);
cMap = zeros(nTrodes,3);
for tr = 1:nTrodes
    [~,tIdx] = min(abs(tMap-repmat(tmpT(tr),length(tMap),1)));
    cMap(tr,:) = cMapInit(tIdx,:);
end


%% which vertices in the MNI surface are closest to electrodes.
% finding distances from each hemisphere.
Xs = vertcat(trodeData.MNItrodeLocs);
for el = 1:nTrodes
    [minDist(el),Vinds(el)] = min(sqrt((Xs(el,1)-vSurf(:,1)).^2+(Xs(el,2)-vSurf(:,2)).^2+(Xs(el,3)-vSurf(:,3)).^2));
end

% plotting a histogram of electrode minimum distances.
plotHistogram = false;
if plotHistogram
    figure(1123)
    subplot(2,1,1)
    hold on
    histogram(minDist, 50,'facecolor','k')
    hold off
    xlabel('distance btw trode and surface (mm)')
    ylabel('count')
    axis square
    
    xLims = xlim;
end

% distance threshold.
% TODO:: compare this to using channels with NMM atlas labels. [20190807]
% TODO:: could also determine this threshold in a data driven way from the
% above histogram [20190809]
distThresh = 9;
% discarding distant electrodes.
elInds = minDist<distThresh;

fprintf('%d total sEEG electrodes.\n%d close enough to brain surface.\nexcluded %d (white matter and extracranial contacts).\n',nTrodes,sum(elInds),nTrodes-sum(elInds))

% plotting effect sizes (t-values) of the linear regression models.
plotEffectSizeScatter = false;
if plotEffectSizeScatter

    figure(1123)
    subplot(2,1,2)
    hold on
    scatter(minDist,tmpT,5,[0 0 0],'filled')
    scatter(minDist(elInds & horzcat(trodeData.pVals)<0.05),tmpT(elInds & horzcat(trodeData.pVals)<0.05),15,cMap(elInds & horzcat(trodeData.pVals)<0.05,:),'filled')
    hold off
    axis square 
    xlim(xLims)

    % saving
    halfMaximize(1123,'left')
    saveas(1123,'~/Dropbox/electrodeDistancesBART_MNI.pdf')
end

% visualization parameters
trodeSize = 20;


%% plotting the contacts for each patient onto the flattened MNI surface.
% if strcmp(whichVis,'2D')
    % now translating projected electrodes from each hemisphere to match
    sEEGverts = vFlat(Vinds(elInds),:);
    modelVerts = vFlat(Vinds(elInds & horzcat(trodeData.pVals)<0.05),:);
    %     sEEGverts(sEEGverts(:,1)<0,1) = sEEGverts(sEEGverts(:,1)<0,1);
    %     sEEGverts(sEEGverts(:,1)>0,1) = sEEGverts(sEEGverts(:,1)>0,1);
    
    pColors = cMap(elInds & horzcat(trodeData.pVals)<0.05,:);
    
    % plotting
%     figure(fHandle)
%     hold on
%     scatter(sEEGverts(:,1),sEEGverts(:,2),round(trodeSize./4),[0 0 0],'filled')
%     scatter(modelVerts(:,1),modelVerts(:,2),trodeSize,pColors,'filled')
%     hold off
% elseif strcmp(whichVis,'3D')
%     %% plotting trodes
%     figure(fHandle)
%     hold on
%     scatter3(ElecXYZMNIRaw(:,1),ElecXYZMNIRaw(:,2),ElecXYZMNIRaw(:,3),trodeSize,[0 1 0],'filled')
%     scatter3(vSurf(Vinds(elInds & [trodeLMM.model_p]<0.05),1),vSurf(Vinds(elInds & [trodeLMM.model_p]<0.05),2),vSurf(Vinds(elInds & [trodeLMM.model_p]<0.05),3),trodeSize,rgb('lightseagreen'),'filled')
%     hold off
%     axis off
% end


%% written summary
%% TODO:: add summary stats for lateralization.
numTrodesTot = length(sEEGverts);
numSigPos = length(vFlat(Vinds(elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)>0),:));
numSigNeg = length(vFlat(Vinds(elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)<0),:));

% per hemisphere results. 
leftHemiIdcs = Xs(:,1)<0;
righHemiIdcs = Xs(:,1)>0;
numSigPosLeft = length(vFlat(Vinds(leftHemiIdcs' & elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)>0),:));
numSigPosRight = length(vFlat(Vinds(righHemiIdcs' & elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)>0),:));
numSigNegLeft = length(vFlat(Vinds(leftHemiIdcs' & elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)<0),:));
numSigNegRight = length(vFlat(Vinds(righHemiIdcs' & elInds & horzcat(trodeData.pVals)<0.05 & horzcat(trodeData.tVals)<0),:));

fprintf('\nTotal number of electrodes: %d, significant Positive effect: %d (%d left, %d right), significnat negative effect: %d (%d left, %d right)',numTrodesTot,numSigPos,numSigPosLeft,numSigPosRight,numSigNeg,numSigNegLeft,numSigNegRight)


%% chi-square Tests 
Xp = [numSigPosLeft numSigPosRight];
Xn = [numSigNegLeft numSigNegRight];
N = [sum(leftHemiIdcs) sum(righHemiIdcs)];

fprintf('\npositive effect lateralization results:\n')
[Hp,pp,chi2p,dfp] = prop_test(Xp,N,false)

fprintf('\nnegative effect lateralization results:\n')
[Hn,pn,chi2n,dfn] = prop_test(Xn,N,false)




