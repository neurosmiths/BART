function [bhvStruct] = BARTbehavior(ptID,plotFlag)
% BARTBEHAVIOR analyzes BART behavior and outputs a detailed struct.
%
% second arg is whether to plot


% ::blackrock data "channel" codes::
%1: trial start ::      [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]
%2: responded ::        [22]
%3: inflating ::        [23 24] = [start stop]
%4: banked ::           [25]
%5: popped ::           [26]
%6: outcome shown ::    [100 101] = [correct incorrect]
%7: max rt exceeded  :: [127]
%8: trial over ::		[120]

% author: EHS20181005
%Edited: RC20220215

% patient details.
ptID = ['202407']; % demo data from one participant
plotFlag = true; % create figures

% loading behavioral data mat file
[BehaviorFile, BehaviorLocation] = uigetfile('.mat', "BehaviorDemoData.mat") % Select Behavior file from repository
matFile =  fullfile(BehaviorLocation,BehaviorFile) % This will allow you to get the "BehaviorDemoData.mat" and save it as matFile.

load(matFile)

% finding nev data to get behavioral markers from neural event file.
[NeuralFile, NeuralLocation] = uigetfile('.nev', "NeuralEventDemoData.mat") % Select Neural file from repository
nevFile =  fullfile(NeuralLocation,NeuralFile) % This will allow you to get the "NeuralEventDemoData.mat" and save it as nevFile.
% trodeLabels = ptTrodesBART(ptID);

% load and define triggers from nevFile
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% [20170713] I made a small error in the handle_input2.m script in the task
% code, such that sometimes there is an additional 23 after the initial 23
% (signifying the start of the inflation). The next two lines remove that
% second 23.
infIdx = trigs==23;
infIdx([diff(infIdx)==0; false] & infIdx==1) = 0;

% response times
respTimes = trigTimes(trigs==25 | trigs==26);
trialStarts = trigTimes(trigs<15);
trialStarts = trialStarts(2:length(respTimes)+1);
% inflateStarts = trigTimes(infIdx);
nTrials = length(respTimes);

% sanity check to compare RTs from blackrock and
sanityCheck = false;
if sanityCheck
    % calculating RTs.
    brRTs = respTimes-trialStarts;
    mlRTs = [data(1:nTrials).inflate_time];
    ctrlIdcs = ~logical([data(1:nTrials).is_control]);
    
    maxRTerror = max(abs(brRTs(ctrlIdcs)-mlRTs(ctrlIdcs)));

    figure(678)
    hold on
    scatter(brRTs(ctrlIdcs),mlRTs(ctrlIdcs))
    plot(1:7,1:7,'--k')
    hold off
    xlabel('RTs calculated from blackrock')
    ylabel('RTs in mat data')
    pause(30)
    close(678)
end


%% [20181009] first organize data
% overall parameters
figID = 1124;
histoAlpha = 0.3;
% in case we want to specify the number of histogram bins.
% uses these many bins for K-L Divergence impulsivity metric.
% leave empty to use the defaults from matlab. ~10-12 is reasonable.
nBins = 10;

% caluclating results.
data = data(1:nTrials);
isCtrl = logical([data.is_control]);
nCtrls = sum(isCtrl);
trialTypeWithCtrls = [data.trial_type];
% 3 = red; 2 = orange; 1 = yellow; 4 = control
greyCtrls = isCtrl & [data.trial_type] == 4; % no points

% 1.1) accuracy
banked = strcmp({data.result},'banked');
popped = strcmp({data.result},'popped');
pTrialBanked = [false, banked(1:end-1)];
pTrialPopped = [false, popped(1:end-1)];
nPopped = sum(popped);

% 1.2) rewarded and unrewarded variables
rewarded = banked & ~greyCtrls; % banked, but take out grey controls.
unrewarded = popped | greyCtrls; % all popped and all grey control trials.
% 1.3) inflation times for rewarded and unrewarded trials
%inflationTimeReward = inflationTime(rewarded);
%inflationTimeUnreward = inflationTime(unrewarded);

% 1.2) Color Accuracy
redBanked = banked & [data.trial_type] == 3;
orangeBanked = banked & [data.trial_type] == 2;
yellowBanked = banked & [data.trial_type] == 1;
redPopped = popped & [data.trial_type] == 3;
orangePopped = popped & [data.trial_type] == 2;
yellowPopped = popped & [data.trial_type] == 1;

% overall accuracy percentages
redAccuracy = (sum(banked(trialTypeWithCtrls==3))./sum(trialTypeWithCtrls==3))*100;
orangeAccuracy = (sum(banked(trialTypeWithCtrls==2))./sum(trialTypeWithCtrls==2))*100;
yellowAccuracy =  (sum(banked(trialTypeWithCtrls==1))./sum(trialTypeWithCtrls==1))*100;


% 2) RT
allRTs = [data.rt];
% Removing Outliers
outRTs = outliers(allRTs); %NIH definition
freeRTs = allRTs(~isCtrl & ~outRTs);
ctrlRTs = allRTs(isCtrl  & ~outRTs);
meanFreeRTs = mean(freeRTs);
bankedRTs = allRTs(banked & ~isCtrl & ~outRTs); % Response times for successful trials (banked)
poppedRTs = allRTs(popped & ~isCtrl & ~outRTs); % Response times for successful trials (popped)
redRTs = allRTs(trialTypeWithCtrls==3 & ~outRTs);
meanRedRTs = mean(redRTs);
orangeRTs = allRTs(trialTypeWithCtrls==2 & ~outRTs);
meanOrangeRTs = mean(orangeRTs);
yellowRTs = allRTs(trialTypeWithCtrls==1 & ~outRTs);
meanYellowRTs = mean(yellowRTs);
redbankedRTs = allRTs(redBanked & ~isCtrl & ~outRTs); % Inflation times for successful RED trials (banked)
redpoppedRTs = allRTs(redPopped & ~isCtrl & ~outRTs); % Inflation times for unsuccessful RED trials (popped)
orangebankedRTs = allRTs(orangeBanked & ~isCtrl & ~outRTs); % Inflation times for successful ORANGE trials (banked)
orangepoppedRTs = allRTs(orangePopped & ~isCtrl & ~outRTs); % Inflation times for unsuccessful ORANGE trials (popped)
yellowbankedRTs = allRTs(yellowBanked & ~isCtrl & ~outRTs); % Inflation times for successful YELLOW trials (banked)
yellowpoppedRTs = allRTs(yellowPopped & ~isCtrl & ~outRTs); % Inflation times for unsuccessful YELLOW trials (popped)
% [20181009] RT isn't really that interesting across conditions. maybe
% just report as text...

postBankRTs = allRTs(~outRTs & pTrialBanked);
postPopRTs = allRTs(~outRTs & pTrialPopped);

pTrialRedBanked = [false, redBanked(1:end-1)];
pTrialRedPopped = [false, redPopped(1:end-1)];
pTrialOrangeBanked = [false, orangeBanked(1:end-1)];
pTrialOrangePopped = [false, orangePopped(1:end-1)];
pTrialYellowBanked = [false, yellowBanked(1:end-1)];
pTrialYellowPopped = [false, yellowPopped(1:end-1)];

postBankRedRTs = allRTs(~outRTs & pTrialRedBanked);
postPopRedRTs = allRTs(~outRTs & pTrialRedPopped);
postBankOrgRTs = allRTs(~outRTs & pTrialOrangeBanked);
postPopOrgRTs = allRTs(~outRTs & pTrialOrangePopped);
postBankYelRTs = allRTs(~outRTs & pTrialYellowBanked);
postPopYelRTs = allRTs(~outRTs & pTrialYellowPopped);



% % save RTs
bhvStruct.allRTs = allRTs;
bhvStruct.freeRTs = freeRTs;
bhvStruct.meanFreeRTs = meanFreeRTs;
bhvStruct.meanRedRTs = meanRedRTs;
bhvStruct.meanOrangeRTs = meanOrangeRTs;
bhvStruct.meanYellowRTs = meanYellowRTs;
bhvStruct.ctrlRTs = ctrlRTs;
bhvStruct.bankedRTs = bankedRTs;
bhvStruct.poppedRTs = poppedRTs;
bhvStruct.redRTs = redRTs;
bhvStruct.orangeRTs = orangeRTs;
bhvStruct.yellowRTs = yellowRTs;
bhvStruct.postPopRTs = postPopRTs;
bhvStruct.postBankRTs = postBankRTs;
% bhvStruct.redbankedRTs = redbankedRTs;
% bhvStruct.orangebankedRTs = orangebankedRTs;
% bhvStruct.yellowbankedRTs = yellowbankedRTs;
% bhvStruct.redpoppedRTs = redpoppedRTs;
% bhvStruct.orangepoppedRTs = orangepoppedRTs;
% bhvStruct.yellowpoppedRTs = yellowpoppedRTs;
% bhvStruct.postBankRedRTs = postBankRedRTs;
% bhvStruct.postPopRedRTs = postPopRedRTs;
% bhvStruct.postBankOrgRTs = postBankOrgRTs;
% bhvStruct.postPopOrgRTs = postPopOrgRTs;
% bhvStruct.postBankYelRTs = postBankYelRTs;
% bhvStruct.postPopYelRTs = postPopYelRTs;
bhvStruct.rewardedTrials = rewarded;

% 3) IT
allITs = [data.inflate_time];
freeITs = allITs(~isCtrl);
ctrlITs = allITs(isCtrl); % both passive and grey control trials
passiveCtrlITs = allITs(trialTypeWithCtrls==1 & isCtrl | trialTypeWithCtrls==2 & isCtrl |trialTypeWithCtrls==3 & isCtrl); % Inflation times for passive balloons only.
greyCtrlITs = allITs(trialTypeWithCtrls==4 & isCtrl); % Inflations for Grey Control Trials only
bankedITs = allITs(banked & ~isCtrl); % Inflation times for successful trials (banked)
poppedITs = allITs(popped & ~isCtrl); % Inflation times for unsuccessful trials (popped)
redITs = allITs(trialTypeWithCtrls==3 & ~isCtrl);
orangeITs = allITs(trialTypeWithCtrls==2 & ~isCtrl);
yellowITs = allITs(trialTypeWithCtrls==1 & ~isCtrl);
redbankedITs = allITs(redBanked & ~isCtrl); % Inflation times for successful RED trials (banked)
redpoppedITs = allITs(redPopped & ~isCtrl); % Inflation times for unsuccessful RED trials (popped)
orangebankedITs = allITs(orangeBanked & ~isCtrl); % Inflation times for successful ORANGE trials (banked)
orangepoppedITs = allITs(orangePopped & ~isCtrl); % Inflation times for unsuccessful ORANGE trials (popped)
yellowbankedITs = allITs(yellowBanked & ~isCtrl); % Inflation times for successful YELLOW trials (banked)
yellowpoppedITs = allITs(yellowPopped & ~isCtrl); % Inflation times for unsuccessful YELLOW trials (popped)
%bhvStruct.inflationTimeReward = inflationTimeReward;
%bhvStruct.inflationTimeUnreward = inflationTimeUnreward;

% Create a normalized IT
normITs = (allITs-mean(allITs))/std(allITs);

postBankITs = normITs(pTrialBanked);
postPopITs = normITs(pTrialPopped);

bhvStruct.poppedTrials = popped; % to look at popped trials over time.

freeTrials = length(freeITs);

% first and last 50 trials ITs for colors.
meanFirstITs = mean((freeITs(1:50)));
meanLastITs = mean((freeITs(freeTrials-49:freeTrials)));
stdFirstITs = std((freeITs(1:50)));
stdLastITs = std((freeITs(freeTrials-49:freeTrials)));

%minumum ITs average for each participant
%minFirstITs = min(freeITs(1:50));
%minLastITs =  min((freeITs(freeTrials-49:freeTrials)));

bhvStruct.meanFirstITs = meanFirstITs;
bhvStruct.meanLastITs = meanLastITs;
bhvStruct.stdFirstITs = stdFirstITs;
bhvStruct.stdLastITs = stdLastITs;

figure(10101234)
hold on
bar(freeITs(1:50),"blue")
bar(freeITs(freeTrials-49:freeTrials),"red")
legend(["First 50 Trials"; "Last 50 Trials"])
text(20,16,sprintf('Mean IT First 50: %s', meanFirstITs));
text(20,15,sprintf('Mean IT Last 50: %s', meanLastITs));
hold off
close(10101234)

% 4) points (proxy for IT)
allPoints = [data.points];
bhvStruct.allPoints = sum(allPoints);
bhvStruct.freePoints = sum(allPoints(~isCtrl));
bhvStruct.ctrlPoints = sum(allPoints(isCtrl));
bhvStruct.redPoints = allPoints(trialTypeWithCtrls==3);
bhvStruct.orangePoints = allPoints(trialTypeWithCtrls==2);
bhvStruct.yellowPoints = allPoints(trialTypeWithCtrls==1);
bhvStruct.freeITs = mean(freeITs); % ADDED IN 02/14
bhvStruct.ctrlITs = mean(ctrlITs);
bhvStruct.bankedITs = bankedITs; % ADDED IN 03/13
bhvStruct.poppedITs = poppedITs;
bhvStruct.redITs = redITs;
bhvStruct.orangeITs = orangeITs;
bhvStruct.yellowITs = yellowITs;
bhvStruct.redbankedITs = redbankedITs;
bhvStruct.orangebankedITs = orangebankedITs;
bhvStruct.yellowbankedITs = yellowbankedITs;
bhvStruct.redpoppedITs = redpoppedITs;
bhvStruct.orangepoppedITs = orangepoppedITs;
bhvStruct.yellowpoppedITs = yellowpoppedITs;
bhvStruct.postBankITs = postBankITs;
bhvStruct.postPopITs = postPopITs;
bhvStruct.passiveCtrlITs = passiveCtrlITs;
bhvStruct.greyCtrlITs = greyCtrlITs;
% bhvStruct.postBankRedITs = postBankRedITs;
% bhvStruct.postPopRedITs = postPopRedITs;
% bhvStruct.postBankOrgITs = postBankOrgITs;
% bhvStruct.postPopOrgITs = postPopOrgITs;
% bhvStruct.postBankYelITs = postBankYelITs;
% bhvStruct.postPopYelITs = postPopYelITs;

% 5) Post-Outcome Slowing (adaptive/maladaptive)
% Need to adjust this for the specific balloon inflation (color) by the
% mean of the rest of the specific balloon colors.

% 5.1) Post-Outcome Slowing (Inflation Times)
% Pop - Bank
%deltaITs = mean(poppedITs) - mean(bankedITs);
% deltaRedITs = mean(redpoppedITs) - mean(redbankedITs);
% deltaOrangeITs =  mean(orangepoppedITs) - mean(orangebankedITs);
% deltaYellowITs =  mean(yellowpoppedITs) - mean(yellowbankedITs);
% % PostPop - PostBank
deltaPostITs = mean(postPopITs) - mean(postBankITs);
% deltaPostRedITs = mean(postPopRedITs) - mean(postBankRedITs);
% deltaPostOrangeITs = mean(postPopOrgITs) - mean(postBankOrgITs);
% deltaPostYellowITs = mean(postPopYelITs) - mean(postBankYelITs); % Change in inflation times. Positive delta indicates an increase in postPop IT.

% 5.2) Post-Outcome Slowing (Response Times)
% Pop - Bank
 deltaRTs = mean(poppedRTs) - mean(bankedRTs);
 deltaRedRTs = mean(redpoppedRTs) - mean(redbankedRTs);
 deltaOrangeRTs =  mean(orangepoppedRTs) - mean(orangebankedRTs);
 deltaYellowRTs =  mean(yellowpoppedRTs) - mean(yellowbankedRTs);
% % PostPop - PostBank
 deltaPostRTs = mean(postPopRTs) - mean(postBankRTs);
 deltaPostRedRTs = mean(postPopRedRTs) - mean(postBankRedRTs);
 deltaPostOrangeRTs = mean(postPopOrgRTs) - mean(postBankOrgRTs);
 deltaPostYellowRTs = mean(postPopYelRTs) - mean(postBankYelRTs); % Change in inflation times. Positive delta indicates an increase in postPop RT.
%
% Calculating the difference between the IT/RT in correct trialsn that either followed an error or a correct choice in trials.
PopToBankITs = normITs(banked &  pTrialPopped); % using norms because ITs vary too much
BankToBankITs = normITs(banked &  pTrialBanked);
PopToBankRTs = allRTs(banked &  pTrialPopped);
BankToBankRTs = allRTs(banked &  pTrialBanked);

bhvStruct.PopToBankITs = PopToBankITs;
bhvStruct.BankToBankITs = BankToBankITs;
bhvStruct.PopToBankRTs = PopToBankRTs;
bhvStruct.BankToBankRTs = BankToBankRTs;

deltaITs = mean(PopToBankITs) - mean(BankToBankITs);%adaptive and maladaptive normITs
deltaRTs = mean(PopToBankRTs) - mean(BankToBankRTs);

bhvStruct.deltaITs = deltaITs;
bhvStruct.deltaPostITs = deltaPostITs;
bhvStruct.deltaRTs = deltaRTs;
bhvStruct.deltaPostRTs = deltaPostRTs;
bhvStruct.deltaRedRTs = deltaRedRTs;
bhvStruct.deltaOrangeRTs = deltaOrangeRTs;
bhvStruct.deltaYellowRTs = deltaYellowRTs;
bhvStruct.deltaPostRedRTs = deltaPostRedRTs;
bhvStruct.deltaPostOrangeRTs = deltaPostOrangeRTs;
bhvStruct.deltaPostYellowRTs = deltaPostYellowRTs;

[P_ITdm,H_ITdm, Stat_ITdm] = ranksum(PopToBankITs, BankToBankITs);
[P_RTdm,H_RTdm, Stat_RTdm] = ranksum(PopToBankRTs, BankToBankRTs);

bhvStruct.inflationDM = H_ITdm;
bhvStruct.responseDM = H_RTdm;

% Simpler version
postRTs = mean(postPopRTs) - mean(postBankRTs);
postITs = mean(postPopITs) - mean(postBankITs);

bhvStruct.postRTs = postRTs;
bhvStruct.postITs = postITs;

[P_postRTdm,H_postRTdm, Stat_postRTdm] = ranksum(postPopRTs, postBankRTs);
[P_postITdm,H_postITdm, Stat_postITdm] = ranksum(postPopITs, postBankITs);

bhvStruct.postInflationDM = H_postITdm;
bhvStruct.postResponseDM = H_postRTdm;

%% [20190311] adding values to the bhvStruct - more later...
bhvStruct.totalTrials = nTrials;
bhvStruct.nCtrls = nCtrls;
bhvStruct.nActiveTrials = nTrials/nCtrls;
bhvStruct.redAccuracy = redAccuracy;
bhvStruct.orangeAccuracy = orangeAccuracy;
bhvStruct.yellowAccuracy = yellowAccuracy;

%% more stats
% histograms for kldiv
maxIT = max(allITs);
nBns = 10;
[ctrlHist,edges] = histcounts(ctrlITs,0:maxIT./nBns:maxIT);
[freeHist,edges] = histcounts(freeITs,0:maxIT./nBns:maxIT);

% Define the bin edges you want
%edgesColor = 1:0.5:10;

% color histograms
[ctrlHistRED,edges] = histcounts(allITs(trialTypeWithCtrls==3 & isCtrl),nBns:maxIT);
[freeHistRED,edges] = histcounts(allITs(trialTypeWithCtrls==3 & ~isCtrl),nBns:maxIT);
[ctrlHistORA,edges] = histcounts(allITs(trialTypeWithCtrls==2 & isCtrl),nBns:maxIT);
[freeHistORA,edges] = histcounts(allITs(trialTypeWithCtrls==2 & ~isCtrl),nBns:maxIT);
[ctrlHistYEL,edges] = histcounts(allITs(trialTypeWithCtrls==1 & isCtrl),nBns:maxIT);
[freeHistYEL,edges] = histcounts(allITs(trialTypeWithCtrls==1 & ~isCtrl),nBns:maxIT);

c = midpoint(edges);
centers = edges(1:end-1)+c(1);
eps = 0.00000000000000000001;
%% ~~~~~~~~~~~ IMPULSIVITY METRIC ~~~~~~~~~~~~~~~~~~~~~~~~~~
% [20220627EHS] more info on this calc:
% distance of active trial distribution from passive trial distribution.
% So lower divergence values means more similar distributions
% => less impulsive choosers.
d=sum((freeHist./sum(freeHist)).*log2(freeHist./sum(freeHist)+eps)-(freeHist./sum(freeHist)).*log2(ctrlHist./sum(ctrlHist)+eps)); %KL(h1,h2)
dSANITY=sum((freeHist./sum(freeHist)).*log2((freeHist./sum(freeHist)+eps)./(ctrlHist./sum(ctrlHist)+eps))); %KL(h1,h2)
bhvStruct.impulsivityKLD = d; % *(sum(banked)./nTrials)
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% [20230214] control for looking at KLD wthin balloon color
dRED=sum((freeHistRED./sum(freeHistRED)).*log2(freeHistRED./sum(freeHistRED)+eps)-(freeHistRED./sum(freeHistRED)).*log2(ctrlHistRED./sum(ctrlHistRED)+eps)); %KL(h1,h2)
dORANGE=sum((freeHistORA./sum(freeHistORA)).*log2(freeHistORA./sum(freeHistORA)+eps)-(freeHistORA./sum(freeHistORA)).*log2(ctrlHistORA./sum(ctrlHistORA)+eps)); %KL(h1,h2)
dYELLOW=sum((freeHistYEL./sum(freeHistYEL)).*log2(freeHistYEL./sum(freeHistYEL)+eps)-(freeHistYEL./sum(freeHistYEL)).*log2(ctrlHistYEL./sum(ctrlHistYEL)+eps)); %KL(h1,h2)
bhvStruct.impulsivityKLD_reds = dRED;
bhvStruct.impulsivityKLD_oranges = dORANGE;
bhvStruct.impulsivityKLD_yellows = dYELLOW;

% saving KLD sanity check stuff: non-parametric comparison of means
[ITp,~,ITstats] = ranksum(ctrlITs,freeITs);
bhvStruct.activeVsPassiveStats.pVal = ITp;
bhvStruct.activeVsPassiveStats.zVal = ITstats.zval;


%% plot a summary figure.
% if the figure exists, close it.
if plotFlag
    if ishandle(figID); close(figID); end
    figure(figID);

    %    text in the top right square
    plotmultipleaxes(1,2,3,0.05,figID);
end

% stats
redBanked = zeros(1,length(banked));
redBanked(trialTypeWithCtrls==3)=1;
redBanked(logical(popped))=0;
orangeBanked = zeros(1,length(banked));
orangeBanked(trialTypeWithCtrls==2)=1;
orangeBanked(logical(popped))=0;
yellowBanked = zeros(1,length(banked));
yellowBanked(trialTypeWithCtrls==1)=1;
yellowBanked(logical(popped))=0;
[accuracyContingencyTbl,chi2Accuracy,pAccuracy] = crosstab(redBanked,orangeBanked,yellowBanked);

if plotFlag
    % plotted text
    title('BART behavior:');
    hold on;
    maxLim = 24;
    text(0,maxLim,sprintf('Subject: %s',ptID));
    text(0,maxLim-2,sprintf('total balloons = %d (%d controls)',nTrials,nCtrls));
    text(0,maxLim-4,sprintf('Accuracy for red balloons: %.1f percent.',redAccuracy));
    text(0,maxLim-6,sprintf('Accuracy for orange balloons: %.1f percent.',orangeAccuracy));
    text(0,maxLim-8,sprintf('Accuracy for yellow balloons: %.1f percent.',yellowAccuracy));
    if pAccuracy<0.01
        text(0,maxLim-10,sprintf('Signifcant difference in accuracy among balloon colors: X^2(%d)=%.1f, p=%.2d',nTrials-2,chi2Accuracy,pAccuracy));
    else
        text(0,maxLim-10,sprintf('NO signifcant difference among balloon colors: X^2(%d)=%.1f, p=%.2d',nTrials-2,chi2Accuracy,pAccuracy));
    end
    text(0,maxLim-12,sprintf('Total accuracy: %.1f percent.',(sum(banked)./nTrials)*100));
    text(0,maxLim-14,sprintf('KLD-derived impulsivity = %.2f ',d));
    text(0,maxLim-16,sprintf('KLD-derived impulsivity (reds only) = %.2f ',dRED));
    text(0,maxLim-18,sprintf('KLD-derived impulsivity (oranges only) = %.2f ',dORANGE));
    text(0,maxLim-20,sprintf('KLD-derived impulsivity (yellows only) = %.2f ',dYELLOW));

    hold off;
    % deets
    ylim([0 25]);
    axis off;
end

%% adding stats to the bhvStruct
bhvStruct.accuracyTot = (sum(banked)./nTrials)*100;
bhvStruct.accuracyStats.df = nTrials-2;
bhvStruct.accuracyStats.chi2Acciracy = chi2Accuracy;
bhvStruct.accuracyStats.pAccuracy = pAccuracy;


%% inflate time histogramsfor free vs. ctrls
% stats
[pColIT,tblColIT,statsColIT] = kruskalwallis(allITs(trialTypeWithCtrls<4 & ~isCtrl),trialTypeWithCtrls(trialTypeWithCtrls<4 & ~isCtrl),'off');

if plotFlag
    plotmultipleaxes(2,2,3,0.05,figID);
    % stats
    title(sprintf('ranksum test, U=%d, p=%.2d.',ITstats.ranksum,ITp));
% plots
hold on
cith = histfit(ctrlITs,nBins);
fith = histfit(freeITs,nBins);

hold off
% deets
cith(1).FaceColor = rgb('forestgreen');
cith(1).FaceAlpha = histoAlpha;
cith(2).Color = rgb('forestgreen');
cith(2).LineWidth = 1;
fith(1).FaceColor = rgb('navy');
fith(1).FaceAlpha = histoAlpha;
fith(2).Color = rgb('navy');
fith(2).LineWidth = 1;
axis square
ylabel('IT count')
xlabel('time (s)')
xlim([0 max(allITs)])

    %% inflate time histograms for balloon types
    plotmultipleaxes(5,2,3,0.05,figID);
    title(sprintf('excluding passive, kruskal-wallis, X^2=%.2f, p=%.2d.',tblColIT{2,5},tblColIT{2,6}));
   % plots
hold on
yith = histfit(yellowITs,nBins);
oith = histfit(orangeITs,nBins);
rith = histfit(redITs,nBins);
hold off
% deets
rith(1).FaceColor = rgb('red');
rith(1).FaceAlpha = histoAlpha;
oith(1).FaceColor = rgb('orangered');
oith(1).FaceAlpha = histoAlpha;
yith(1).FaceColor = rgb('gold');
yith(1).FaceAlpha = histoAlpha;
rith(2).Color = rgb('red');
rith(2).LineWidth = 1;
oith(2).Color = rgb('orangered');
oith(2).LineWidth = 1;
yith(2).Color = rgb('gold');
yith(2).LineWidth = 1;
axis square
ylabel('IT count')
xlabel('time (s)')
xlim([0 max(allITs)])
end

%% more stats
bhvStruct.balloonColorInflationTimeStats.pVal = pColIT;
bhvStruct.balloonColorInflationTimeStats.statTbl = tblColIT;
bhvStruct.balloonColorInflationTimeStats.statStruct = statsColIT;


%% [20181009] plot [RT] distributions as violin plots
% stats
[RTp,~,RTstats] = ranksum(ctrlRTs,freeRTs);

if plotFlag
    plotmultipleaxes(15,7,3,0.05,figID);
    % plots
    hold on
    violinPlot(ctrlRTs','showMM',4,'color',rgb('forestgreen'),'xValues',1)
    violinPlot(freeRTs','showMM',4,'color',rgb('navy'),'xValues',2)
    hold off
    % deets
    ylabel('RT (s)')
    xlabel(sprintf('ranksum, U(2)=%d, p=%.2f.',RTstats.ranksum,RTp))
    if ceil(max(allRTs))>10
        ylim([0 10])
    else
        ylim([0 ceil(max(allRTs))])
    end
end

%% more stats
bhvStruct.reactionTimeStats.pVal = RTp;
bhvStruct.reactionTimeSTats.zVal = RTstats.zval;

%% [20190116] plot [RT] distributions as violin plots
% stats
[pColRT,tblColRT,statsColRT] = kruskalwallis(allRTs(trialTypeWithCtrls<4 & ~isCtrl),trialTypeWithCtrls(trialTypeWithCtrls<4 & ~isCtrl),'off');

if plotFlag
    plotmultipleaxes(18,7,3,0.05,figID);
    % plots
    hold on
    violinPlot(redRTs','showMM',4,'color',rgb('red'),'xValues',1)
    violinPlot(orangeRTs','showMM',4,'color',rgb('darkorange'),'xValues',2)
    violinPlot(yellowRTs','showMM',4,'color',rgb('gold'),'xValues',3)
    hold off
    % deets
    ylabel('RT (s)')
    xlabel(sprintf('kruskal-wallis, X^2(2)=%.2f, p=%.2f.',tblColRT{2,5},tblColRT{2,6}))
    if ceil(max(allRTs))>10
        ylim([0 10])
    else
        ylim([0 ceil(max(allRTs))])
    end
end

%% more stats
bhvStruct.balloonColorRTStats.pVal = pColRT;
bhvStruct.balloonColorRTStats.statStruct = statsColRT;

%% [20181009] plot follownig trial [RT] distributions as violin plots
% stats
[RTp,~,RTstats] = ranksum(postBankRTs,postPopRTs);

if plotFlag
    plotmultipleaxes(21,7,3,0.05,figID);
    % plots
    hold on
    violinPlot(postBankRTs','showMM',4,'color',rgb('forestgreen'),'xValues',1)
    violinPlot(postPopRTs','showMM',4,'color',rgb('darkorange'),'xValues',2)
    hold off
    % deets
    ylabel('RT (s)')
    xlabel(sprintf('ranksum, U(2)=%d, p=%.2f.',RTstats.ranksum,RTp))
    if ceil(max(allRTs))>10
        ylim([0 10])
    else
        ylim([0 ceil(max(allRTs))])
    end
end

%% more stats
bhvStruct.postBankPopRTStats.pVal = RTp;
bhvStruct.postBankPopRTStats.zVal = RTstats.zval;
bhvStruct.postBankPopRTStats.bankPopMeans = [mean(postBankRTs) mean(postPopRTs)];
bhvStruct.postBankPopRTStats.bankPopSTDs = [std(postBankRTs) std(postPopRTs)];

% plotting points histograms for balloon types
% stats
[pColPts,tblColPts,statsColPTs] = kruskalwallis(allPoints(trialTypeWithCtrls<4),trialTypeWithCtrls(trialTypeWithCtrls<4),'off');

if plotFlag
    plotmultipleaxes(3,2,3,0.05,figID);
    title(sprintf('including passive, kruskal-wallis, X^2=%.2f, p=%.2d.',tblColPts{2,5},tblColPts{2,6}))
    % plots
    hold on
    yith = histfit(bhvStruct.yellowPoints,nBins, 'gamma');
    oith = histfit(bhvStruct.orangePoints,nBins, 'gamma');
    rith = histfit(bhvStruct.redPoints,nBins, 'gamma');
    hold off

    % deets
    rith(1).FaceColor = rgb('red');
    rith(1).FaceAlpha = histoAlpha;
    oith(1).FaceColor = rgb('darkorange');
    oith(1).FaceAlpha = histoAlpha;
    yith(1).FaceColor = rgb('gold');
    yith(1).FaceAlpha = histoAlpha;
    rith(2).Color = rgb('red');
    rith(2).LineWidth = 1;
    oith(2).Color = rgb('darkorange');
    oith(2).LineWidth = 1;
    yith(2).Color = rgb('gold');
    yith(2).LineWidth = 1;
    axis square
    ylabel('point bin count')
    xlabel('points')
    xlim([0 max(allPoints)])

end

%% More behavioral modeling for each patient.
% - model current inflation time as a fucntion of 1) previous trial
% inflation time,
% build table

% trialTypeCategorical = categorical(trialTypeWithCtrls);
% trialTypeCategorical(trialTypeWithCtrls==1)= 'A';
% trialTypeCategorical(trialTypeWithCtrls==2)= 'B';
% trialTypeCategorical(trialTypeWithCtrls==3)= 'C';
% trialTypeCategorical(trialTypeWithCtrls==4)= 'D';

pointsTable_wCTRL = table([data.points]',[data(2:end).points NaN]',trialTypeWithCtrls',...
    'VariableNames',{'points','previousPoints','balloonType'});
pointsTable = pointsTable_wCTRL(~isCtrl,:);

%bhvStruct.glmPP = fitglme(pointsTable,'points ~ 1 + previousPoints + (1 + previousPoints|balloonType)') %  + (1 + 1|balloonType)
% bhvStruct.glmPP_p = bhvStruct.glmPP.Coefficients{2,6};
% bhvStruct.glmPP_t = bhvStruct.glmPP.Coefficients{2,4};
% bhvStruct.glmPP_R2 =  bhvStruct.glmPP.Rsquared.ordinary;

% bhvStruct.glmBT = fitglme(pointsTable,'points ~ 1 + balloonType + (1 + 1|balloonType)')
% bhvStruct.glmBT_p = bhvStruct.glmBT.Coefficients{2,6};
% bhvStruct.glmBT_t = bhvStruct.glmBT.Coefficients{2,4};
% bhvStruct.glmBT_R2 =  bhvStruct.glmBT.Rsquared.ordinary;
%
% bhvStruct.LMMcomparison = compare(bhvStruct.glmPP,bhvStruct.glmBT)

%% RC: More behavioral modeling for each patient.
% - model current inflation time as a function of 1) same trial
% response time,
% build table

% trialTypeCategorical = categorical(trialTypeWithCtrls);
% trialTypeCategorical(trialTypeWithCtrls==1)= 'A';
% trialTypeCategorical(trialTypeWithCtrls==2)= 'B';
% trialTypeCategorical(trialTypeWithCtrls==3)= 'C';
% trialTypeCategorical(trialTypeWithCtrls==4)= 'D';
%----------
%pointsTable_wCTRL = table([data.points]',[data(2:end).points NaN]',trialTypeWithCtrls',...
%    'VariableNames',{'points','previousPoints','balloonType', 'freeRT'});
%pointsTable = pointsTable_wCTRL(~isCtrl,:);

%bhvStruct.glmPP = fitglme(pointsTable,'points ~ 1 + previousPoints + (1 + previousPoints|balloonType)') %  + (1 + 1|balloonType)

% RC: bhvStruct.glmPP = fitglme(pointsTable,'points ~ 1 + freeRT + (1 + freeRT|balloonType)') + (1 + 1|balloonType)
%----------
% bhvStruct.glmPP_p = bhvStruct.glmPP.Coefficients{2,6};
% bhvStruct.glmPP_t = bhvStruct.glmPP.Coefficients{2,4};
% bhvStruct.glmPP_R2 =  bhvStruct.glmPP.Rsquared.ordinary;

% bhvStruct.glmBT = fitglme(pointsTable,'points ~ 1 + balloonType + (1 + 1|balloonType)')
% bhvStruct.glmBT_p = bhvStruct.glmBT.Coefficients{2,6};
% bhvStruct.glmBT_t = bhvStruct.glmBT.Coefficients{2,4};
% bhvStruct.glmBT_R2 =  bhvStruct.glmBT.Rsquared.ordinary;
%
% bhvStruct.LMMcomparison = compare(bhvStruct.glmPP,bhvStruct.glmBT)

% Adding in sex and age data
%try
%SexAge = readtable('D:\Data\preProcessed\BART_preprocessed\BART_Subject_Sex_Age.csv','format','%s%s%D%s');
%SexAgeIdx = contains(SexAge.patientID, ptID);
%bhvStruct.Sex = SexAge.Sex(SexAgeIdx);
%bhvStruct.DOB = SexAge.DOB(SexAgeIdx);
%bhvStruct.Age = SexAge.Age(SexAgeIdx);
%end
