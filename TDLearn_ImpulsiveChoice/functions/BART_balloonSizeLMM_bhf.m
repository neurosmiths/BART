function [trodeLMM] = BART_balloonSizeLMM_bhf(ptID,BARTdir)
% BART_COLORTRACKING_BHF fit linear model to BHF during balloon inflation
%
%   [BARTstats] = BART_colortracking_bhf(ptID,nevFile) analyzes LFP data
%   for the patient specified in the string ptID using the data in nevFile.
%

% author: EHS20191307

ptID = '201911'; % this line only for debugging

if nargin<2
    BARTdir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data'];
end
nevList = dir([BARTdir '/*.nev']);

% loading PTB BHV
load(fullfile(BARTdir,[ptID '.bartBHV.mat']))

% electrode information.
[trodeLabels,isECoG] = ptTrodesBART(ptID);

% colormap
cMap(1,:) = [1 0.9 0];
cMap(2,:) = [1 0.5 0];
cMap(3,:) = [1 0 0];
cMap(4,:) = [0.5 0.5 0.5];

% load and define triggers from nevFle
if length(nevList)>2
    error('2manyNEVs!')
else
    nevFile = fullfile(BARTdir,nevList(1).name);
end
NEV = openNEV(fullfile(BARTdir,nevList(1).name),'overwrite');
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
pre = 1;
post = 1;

% % task parameters (out of order)
% alignName = 'banksandpops';
% bankTimes = trigTimes(trigs==25);
% popTimes = trigTimes(trigs==26);
% nBanks = length(bankTimes);
% nPops = length(popTimes);
% outcomeType = [ones(1,nBanks) 2*ones(1,nPops)]; % 1 = bank, 2 = pop
% nTrials = length(outcomeType);
% %1: trial start ::      [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]

% task parameters in chronological order..
balloonTimes = trigTimes(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
inflateTimes = trigTimes(trigs==23);
balloonType = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14); % 1 = bank, 2 = pop
respTimes = trigTimes(trigs==26 | trigs==25);

% task identifiers
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==4 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
isCTRL = balloonIDs>10;

% only including complete trials.
nTrials = length(respTimes);
balloonType = balloonType(1:nTrials);
inflateTimes = inflateTimes(1:nTrials);
balloonTimes = balloonTimes(1:nTrials);
isCTRL = isCTRL(1:nTrials);

% BHF filter
[b,a] = butter(4,[70 150]/(Fs/2));

if isequal(ptID,'201903')
    [b120,a120] = iirnotch(120/(Fs/2),60/(Fs/2)/35);
end

whichModel = 'passive_each';
% for 'each' electrode or 'all' electrodes,
% or 'active_each' looks at the interaction of active and passive trials.


%% [20196007] new strategy:
% get the BHF signal between inflation start and pop.
% examine trial-by-trial covariation with baloon size.
% linear model.
% - adding a second of data on each edge for safety.

% starting by collecting the data.
for tt = 1:nTrials
    % adding a second of data on each side
    if isequal(ptID,'201903')
        % notch filtering @ 120 Hz
        eachBalloonLFP{tt} = filtfilt(b120,a120,double(NSX.Data(1:nChans,floor(Fs*inflateTimes(tt))-Fs*pre:ceil(Fs*respTimes(tt))+Fs*post)));
    else
        eachBalloonLFP{tt} = NSX.Data(1:nChans,floor(Fs*inflateTimes(tt))-Fs*pre:ceil(Fs*respTimes(tt))+Fs*post);
    end
    
    for ch = 1:nChans
        tmpHG(ch,:) = abs(hilbert(filtfilt(b,a,double(eachBalloonLFP{tt}(ch,:)))));
    end
    
    % now removin the second of data on each side.
    HGcell{tt} = tmpHG(:,Fs*pre:end-Fs*post);
    baloonSizeSurrogate{tt} = 1:size(HGcell{tt},2);
    baloonSizeSurrogate_01{tt} = linspace(0,1,size(HGcell{tt},2));
    trialSurrogate{tt} = tt.*ones(1,size(HGcell{tt},2));
    if strcmp(whichModel,'passive_each')
        isCTRLsurro{tt} = repmat(~isCTRL(tt),1,size(HGcell{tt},2));
    elseif strcmp(whichModel,'active_each')
        isCTRLsurro{tt} = repmat(isCTRL(tt),1,size(HGcell{tt},2));
    end
    updateUser('stacked and filtered data for trial',tt,50,nTrials)
    clear tmpHG
end

if strcmp(whichModel,'each')
    % balloon size and trial
    normBalloonSizeVec = cell2mat(baloonSizeSurrogate_01)';
    balloonSizeVec = cell2mat(baloonSizeSurrogate)';
    
    trialVec = categorical(cell2mat(trialSurrogate)');
    
    % data matrix
    HGmat = cell2mat(HGcell)';
    
    % electrode label vector
    trodeLabelVec = trodeLabels(isECoG);
    
    % saving model results for all electrodes
    
    % looping over electrodes.
    for el = 1:nChans
        fprintf('\nrunning model for electrode %s',trodeLabelVec{el})
        
        % build table.
        HGtbl = table(balloonSizeVec,normBalloonSizeVec,trialVec,HGmat(:,el),'VariableNames',{'balloonSize','normalizedBalloonSize','trial','BHFpower'});
        
        try
            % run statistical model.
            tic
            fprintf('\nrunning linear model. Started at %s.',char(datetime('now')))
            trodeLMM(el).glme = fitglme(HGtbl,'balloonSize ~ 1 + BHFpower + (1 + BHFpower|trial)');
            trodeLMM(el).glme
            trodeLMM(el).model_df = trodeLMM(el).glme.Coefficients{2,5};
            trodeLMM(el).model_t = trodeLMM(el).glme.Coefficients{2,4};
            trodeLMM(el).model_p = trodeLMM(el).glme.Coefficients{2,6};
            trodeLMM(el).model_pCorr = trodeLMM(el).glme.Coefficients{2,6}*nChans;
            trodeLMM(el).trodeLabel = trodeLabelVec{el};
            A = toc;
            fprintf('\nfitting linear model for %s took %.2f minutes',trodeLabelVec{el},A./60)
        catch
            fprintf('\nmodel fitting failed for %s',trodeLabelVec{el})
            trodeLMM(el).glme = 'fail';
            trodeLMM(el).model_df = 0;
            trodeLMM(el).model_t = 0;
            trodeLMM(el).model_p = 1;
            trodeLMM(el).model_pCorr = 1;
            trodeLMM(el).trodeLabel = trodeLabelVec{el};
        end
    end
    
    if ~exist(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'file')
        % saving data if files don't already exist.
        save(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'trodeLMM','-v7.3')
    end
    
    %  plotting examples here...
    figure
    subplot(1,2,1)
    [~,bigT] = max([trodeLMM.model_t]);
    fixedEffectsBig = designMatrix(trodeLMM(bigT).glme,'fixed');
    scatter(HGtbl.normalizedBalloonSize,fixedEffectsBig(:,2),5,[0 0 0],'filled')
    axis square
    xlabel('normalized balloon size')
    ylabel('GLME fixed effects')
    title(sprintf('largest effect: %s',trodeLabelVec{bigT}))
    
    subplot(1,2,2)
    [~,lilT] = min([trodeLMM.model_t]);
    fixedEffectsLil = designMatrix(trodeLMM(lilT).glme,'fixed');
    scatter(HGtbl.normalizedBalloonSize,fixedEffectsLil(:,2),5,[0 0 0],'filled')
    axis square
    xlabel('normalized balloon size')
    ylabel('GLME fixed effects')
    title(sprintf('smallest effect: %s',trodeLabelVec{lilT}))
    
elseif strcmp(whichModel,'active_each')
    
    % balloon size and trial
    normBalloonSizeVec = cell2mat(baloonSizeSurrogate_01)';
    balloonSizeVec = cell2mat(baloonSizeSurrogate)';
    
    trialVec = categorical(cell2mat(trialSurrogate)');
    ctrlVec = categorical(cell2mat(isCTRLsurro)');
    
    % data matrix
    HGmat = cell2mat(HGcell)';
    
    % electrode label vector
    trodeLabelVec = trodeLabels(isECoG);
    
    % saving model results for all electrodes
    if ~exist(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'file')
        
        % looping over electrodes.
        for el = 1:nChans
            fprintf('\nrunning model for electrode %s',trodeLabelVec{el})
            
            % build table.
            HGtbl = table(balloonSizeVec,normBalloonSizeVec,trialVec,ctrlVec,HGmat(:,el),'VariableNames',{'balloonSize','normalizedBalloonSize','trial','isCTRL','BHFpower'});
            
            try
                % run statistical model.
                tic
                fprintf('\nrunning linear model. Started at %s.',char(datetime('now')))
                trodeLMM(el).glme = fitglme(HGtbl,'balloonSize ~ 1 + isCTRL*BHFpower + (1 + isCTRL*BHFpower|trial)');
                trodeLMM(el).glme
                trodeLMM(el).model_df = trodeLMM(el).glme.Coefficients{2,5};
                trodeLMM(el).model_t = trodeLMM(el).glme.Coefficients{2,4};
                trodeLMM(el).model_p = trodeLMM(el).glme.Coefficients{2,6};
                trodeLMM(el).model_pCorr = trodeLMM(el).glme.Coefficients{2,6}*nChans;
                trodeLMM(el).trodeLabel = trodeLabelVec{el};
                A = toc;
                fprintf('\nfitting linear model for %s took %.2f minutes',trodeLabelVec{el},A./60)
            catch
                fprintf('\nmodel fitting failed for %s',trodeLabelVec{el})
                trodeLMM(el).glme = 'fail';
                trodeLMM(el).model_df = 0;
                trodeLMM(el).model_t = 0;
                trodeLMM(el).model_p = 1;
                trodeLMM(el).model_pCorr = 1;
                trodeLMM(el).trodeLabel = trodeLabelVec{el};
            end
        end
        
        % saving data if files don't already exist.
        save(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'trodeLMM','-v7.3')
    else
        % loading data if it already exists
        load(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)))
    end
    
elseif strcmp(whichModel,'passive_each')
    
    % balloon size and trial
    
    normBalloonSizeVec = cell2mat(baloonSizeSurrogate_01)';
    balloonSizeVec = cell2mat(baloonSizeSurrogate)';
    
    trialVec = categorical(cell2mat(trialSurrogate)');
    ctrlVec = categorical(cell2mat(isCTRLsurro)'); % THIS IS THE ONLY CHANGE FROM THE PREVIOUS CASE: ~
    
    % data matrix
    HGmat = cell2mat(HGcell)';
    
    % electrode label vector
    trodeLabelVec = trodeLabels(isECoG);
    
    % saving model results for all electrodes
    if ~exist(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'file')
        
        % looping over electrodes.
        for el = 1:nChans
            fprintf('\nrunning model for electrode %s',trodeLabelVec{el})
            
            % build table.
            HGtbl = table(balloonSizeVec,normBalloonSizeVec,trialVec,ctrlVec,HGmat(:,el),'VariableNames',{'balloonSize','normalizedBalloonSize','trial','isCTRL','BHFpower'});
            
            try
                % run statistical model.
                tic
                fprintf('\nrunning linear model. Started at %s.',char(datetime('now')))
                trodeLMM(el).glme = fitglme(HGtbl,'balloonSize ~ 1 + isCTRL*BHFpower + (1 + isCTRL*BHFpower|trial)');
                trodeLMM(el).glme
                trodeLMM(el).model_df = trodeLMM(el).glme.Coefficients{2,5};
                trodeLMM(el).model_t = trodeLMM(el).glme.Coefficients{2,4};
                trodeLMM(el).model_p = trodeLMM(el).glme.Coefficients{2,6};
                trodeLMM(el).model_pCorr = trodeLMM(el).glme.Coefficients{2,6}*nChans;
                trodeLMM(el).trodeLabel = trodeLabelVec{el};
                A = toc;
                fprintf('\nfitting linear model for %s took %.2f minutes',trodeLabelVec{el},A./60)
            catch
                fprintf('\nmodel fitting failed for %s',trodeLabelVec{el})
                trodeLMM(el).glme = 'fail';
                trodeLMM(el).model_df = 0;
                trodeLMM(el).model_t = 0;
                trodeLMM(el).model_p = 1;
                trodeLMM(el).model_pCorr = 1;
                trodeLMM(el).trodeLabel = trodeLabelVec{el};
            end
        end
        
        % saving data if files don't already exist.
        save(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)),'trodeLMM','-v7.3')
    else
        % loading data if it already exists
        load(fullfile(BARTdir,sprintf('%s_balloonSizeTrackingGLMresults_%sTrodesModel.mat',ptID,whichModel)))
    end
    
    
    % running an omnibus model.
elseif strcmp(whichModel,'all')
    %% [20190307] making a data table
    % balloon size and trial number vectors
    tic
    balloonSizeVec = repmat(cell2mat(baloonSizeSurrogate)',nChans,1);
    trialVec = categorical(repmat(cell2mat(trialSurrogate)',nChans,1));
    
    % broadband high frequency data vectors
    trodeData = HGmat(:);
    
    % electrode label vector
    trodeLabelMat = repmat(trodeLabels(isECoG),1,length(HGmat))';
    trodeLabelVec = trodeLabelMat(:);
    
    % build table.
    HGtbl = table(balloonSizeVec,trialVec,trodeData,trodeLabelVec,'VariableNames',{'balloonSize','trial','BHFpower','electrode'});
    A = toc;
    fprintf('\nbuilding data table took %.2f minutes',A./60)
    
    % run statistical model.
    tic
    fprintf('\nrunning linear model. Started at %s.',char(datetime('now')))
    hglme = fitglme(HGtbl,'balloonSize ~ 1 + BHFpower + (1 + BHFpower|trial) + (1 + BHFpower|electrode)');
    A = toc;
    fprintf('\nfitting linear model took %.2f minutes',A./60)
    
    save(fullfile(BARTdir,sprintf('%s_colorTrackingGLMresults_allTrodesModel.m',ptID,whichModel)),'HGtbl','hglme','-v7.3')
end


