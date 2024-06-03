function [bhvStruct, bisbasStruct] = bisbasBART()

% author: RLC02152023

BisBasptArray = {'201901', '201902', '201905', '201909', '201910', '201911', '201913', '201914', '202001', '202006','202007','202008','202009','202011','202014','202015','202107',...
    '202110','202114','202117','202118','202201','202207','202209','202212','202214','202215','202216','202217','202302'}

% only need patients up to 202302 for TD paper: '202307', '202308' '202314'
% '202006u' (problems with the u)

fprintf('\nnumber of subjects: %d\n',length(BisBasptArray));

% DO I NEED THIS?
  %  ptID = ['202308'];

    % DO I NEED THIS?
  %  parentDir = ['D:\Data\preProcessed\BART_preprocessed\' ptID '\Data\' ptID '.bartBHV.mat'];

    % loading bisbas mastersheet
    BisBas = readtable('D:\Data\preProcessed\BART_preprocessed\BISBAS_TDproject.csv')

   % BisBas.patientID = BisBasptArray
   % BisBas = convertvars(BisBas,@iscell, 'string')
   % BisBas.patientID = cellfun(@(x) str2double(strrep(x, 'NaN', '202006u')), BisBas.patientID);
   % [BisBas] = readcell('D:\Data\preProcessed\BART_preprocessed\BISBAS.xlsx'); % read in as a cell to include NaN patientIDs.
   % BisBas = cell2table(BisBas) % changing 
    BisBas.patientID = num2str(BisBas.patientID); % Make ptID a string variable
    BisBas.patientID = cellstr(BisBas.patientID); % Make ptID a string cell variable


    %% ----------- BISBAS ITEMS ----------- %%
    % Items other than 2 and 22 are reverse-scored.
    % BAS Drive:  3, 9, 12, 21
    % BAS Fun Seeking:  5, 10, 15, 20
    % BAS Reward Responsiveness:  4, 7, 14, 18, 23
    % BIS:  2, 8, 13, 16, 19, 22, 24
    % Items 1, 6, 11, 17,  are fillers.

    % Reverse-Scoring Items ---
    % Items are scored from 1 - 4. To reverse score minus 5 from score.
    % Do not include items 2 & 22, or fillers.
    BisBas.item_3 = 5 - BisBas.item_3;
    BisBas.item_4 = 5 - BisBas.item_4;
    BisBas.item_5 = 5 - BisBas.item_5;
    BisBas.item_7 = 5 - BisBas.item_7;
    BisBas.item_8 = 5 - BisBas.item_8;
    BisBas.item_9 = 5 - BisBas.item_9;
    BisBas.item_10 = 5 - BisBas.item_10;
    BisBas.item_12 = 5 - BisBas.item_12;
    BisBas.item_13 = 5 - BisBas.item_13;
    BisBas.item_14 = 5 - BisBas.item_14;
    BisBas.item_15 = 5 - BisBas.item_15;
    BisBas.item_16 = 5 - BisBas.item_16;
    BisBas.item_18 = 5 - BisBas.item_18;
    BisBas.item_19 = 5 - BisBas.item_19;
    BisBas.item_20 = 5 - BisBas.item_20;
    BisBas.item_21 = 5 - BisBas.item_21;
    BisBas.item_23 = 5 - BisBas.item_23;
    BisBas.item_24 = 5 - BisBas.item_24;

    % Drop Fillers from table
    BisBas = removevars(BisBas,["item_1", "item_6", "item_11", "item_17"]);

    % TOTAL BISBAS SCORE ---
    BisBasVars = ["item_2", "item_3", "item_4", "item_5", "item_7","item_8",...
        "item_9","item_10", "item_12", "item_13","item_14", "item_15", "item_16",...
        "item_18", "item_19","item_20", "item_21", "item_22", "item_23","item_24"];

    BisBas.TotalScore = sum(BisBas{:,BisBasVars},2);   % Total Scores per pt

    % BAS DRIVE ---
    BasDriveVars = ["item_3", "item_9", "item_12", "item_21"];
    BisBas.BasDriveScore = sum(BisBas{:,BasDriveVars},2);  % Drive Scores per pt

    % BAS FUN SEEKING ---
    BasFunVars = ["item_5", "item_10", "item_15", "item_20"];
    BisBas.BasFunScore = sum(BisBas{:,BasFunVars},2); % Fun Scores per pt

    % BAS REWARD RESPONSIVENESS ---
    BasRewVars = ["item_4", "item_7", "item_14", "item_18", "item_23"];
    BisBas.BasRewScore = sum(BisBas{:,BasRewVars},2); % Reward Scores per pt

    % BIS ---
    BisVars = ["item_2", "item_8", "item_13", "item_16", "item_19", "item_22", "item_24"];
    BisBas.BisScore = sum(BisBas{:,BisVars},2); % BIS Scores per pt

    % Across subjects descriptives ----
    BasDriveMean = mean(BisBas.BasDriveScore);
    BasDriveStd = std(BisBas.BasDriveScore);
    BasFunMean = mean(BisBas.BasFunScore);
    BasFunStd = std(BisBas.BasFunScore);
    BasRewMean = mean(BisBas.BasRewScore);
    BasRewStd = std(BisBas.BasRewScore);
    BisMean = mean(BisBas.BisScore);
    BisStd = std(BisBas.BisScore);

    %% ----------- BISBAS CORRELATES WITH KLD SCORES ----------- %%

    % load in BART bhv Struct.
    [ptArray,bhvStruct,hazEEG] = BARTnumbers;

    % looping over patients
    nPts = length(BisBasptArray);

    % Getting impulsivity scores for BISBAS pts
    % median split in KLD metric
    impulsivityMetric = log10([bhvStruct(:).impulsivityKLD]);

    % Gaussian mixture modeling on KLD metric.
    GM = fitgmdist(impulsivityMetric',2);
    idcs = cluster(GM,impulsivityMetric');
    tmp = logical(idcs-1);

    % [20220627EHS] since the GMM doens't classify in order of means, we
    % have to do that here...
    muA = mean(impulsivityMetric(tmp));
    muB = mean(impulsivityMetric(~tmp));
    if muA>muB
        impulsiveChoosers = tmp;
    else
        impulsiveChoosers = ~tmp;
    end

    medSplit = impulsiveChoosers;

    impulsivityBoundary = (min(impulsivityMetric(impulsiveChoosers))+max(impulsivityMetric(~impulsiveChoosers)))/2;


    % ADDING TO
    bhvTable = struct2table(bhvStruct); % Make bhv a table
    bhvTable.patientID = ptArray(:);
    bhvTable.impulsivityMetric = impulsivityMetric(:); % Add IM to bhvTable
    bhvTable.groupSplit = medSplit(:); % 1 = more impulsive

    %keyboard % when you add new patients you'll need to edit the next few lines of code to remove the right patients from the table.

    %subsetting bhv by BisBasptArray
    BBBhv = bhvTable(contains(bhvTable.patientID, BisBasptArray), :); % Choosing patientIDs from substring (bisbasptarray) that appears in bhvTable)
    BBBhv(11,:) = [] % FOR NOW removing patient 202006u (eventually this will be deleted)
    BBBhv(26,:) = [] % FOR NOW removing patient 202212b (eventually this will be deleted)

    BisBasBhv = join(BBBhv, BisBas) % Combining BHV data and BISBAS items
    bisbasBhvStruct = table2struct(BisBasBhv); % Make bisbas a structure

  
%% PLOTS FOR BISBAS vs. IMPULSIVITY %%

impulsiveChoosers_BisBas = BisBasBhv.groupSplit == 1 % finding impulsiveChoosers for BISBAS subset.
impulsivityMetric_BisBas = BisBasBhv.impulsivityMetric % finding impulsivityMetric for BISBAS subset.

%Indexing for Plots
totalScore = [bisbasBhvStruct.TotalScore];
basD = [bisbasBhvStruct.BasDriveScore];
basF = [bisbasBhvStruct.BasFunScore];
basR = [bisbasBhvStruct.BasRewScore];
bisScore = [bisbasBhvStruct.BisScore];
groupSplit = [bisbasBhvStruct.groupSplit];

% DO STATS HERE (BETWEEN GROUPS)
% anova (complete)
X = [totalScore' basD' basF' basR' bisScore'];
[pA,tabA,statsA] = kruskalwallis(X);
[c,m,h,nms] = multcompare(statsA,'display','on');

% ranksum
[ptotalScore,htotalScore,stattotalScore] = ranksum(totalScore(~impulsiveChoosers_BisBas),totalScore(impulsiveChoosers_BisBas)) % NS
[pbasD,hbasD,statbasD] = ranksum(basD(~impulsiveChoosers_BisBas),basD(impulsiveChoosers_BisBas)) % NS
[pbasF,hbasF,statbasF] = ranksum(basF(~impulsiveChoosers_BisBas),basF(impulsiveChoosers_BisBas)) % NS
[pbasR,hbasR,statbasR] = ranksum(basR(~impulsiveChoosers_BisBas),basR(impulsiveChoosers_BisBas)) % NS
[pbisScore,hbisScore,statbisScore] = ranksum(bisScore(~impulsiveChoosers_BisBas),bisScore(impulsiveChoosers_BisBas)) % NS

%plotting
BOXPLOT = true;
figure(1)
subplot(2,3,[1 4])
hold on
line([0 nPts],[impulsivityBoundary impulsivityBoundary],'linestyle','-','color',rgb('indigo'))
text(nPts-2,median(impulsivityMetric_BisBas)+1,'median split','color',rgb('indigo'))
scatter(find(impulsiveChoosers_BisBas),impulsivityMetric_BisBas(impulsiveChoosers_BisBas),20,[0 0 0],'filled','^')
scatter(find(~impulsiveChoosers_BisBas),impulsivityMetric_BisBas(~impulsiveChoosers_BisBas),20,[0 0 0],'filled')
hold off
xlim([0.9 nPts]) %from .7
ylabel('impulsivity metric')
%xlabel('patient number')
title('BISBAS impulsivity: - low: o  - high: ^')

subplot(2,3,[2 5])
hold on
if BOXPLOT
    betterBoxplot(1,basD(impulsiveChoosers_BisBas),rgb('orangered'),15,'^',1,false)
    betterBoxplot(1.5,basD(~impulsiveChoosers_BisBas),rgb('orangered'),15,'o',1,false) % ARE THESE CORRECT (group split for just bisbasbhv)
    betterBoxplot(2,basF(impulsiveChoosers_BisBas),rgb('red'),15,'^',1,false)
    betterBoxplot(2.5,basF(~impulsiveChoosers_BisBas),rgb('red'),15,'o',1,false)
    betterBoxplot(3,basR(impulsiveChoosers_BisBas),rgb('black'),15,'^',1,false)
    betterBoxplot(3.5,basR(~impulsiveChoosers_BisBas),rgb('black'),15,'o',1,false)
    betterBoxplot(4,bisScore(impulsiveChoosers_BisBas),rgb('green'),15,'^',1,false)
    betterBoxplot(4.5,bisScore(~impulsiveChoosers_BisBas),rgb('green'),15,'o',1,false)

else
    scatter(1*ones(1,sum(impulsiveChoosers_BisBas)),basD(impulsiveChoosers_BisBas),10,rgb('orangered'),'filled','^')
    scatter(1.3*ones(1,sum(~impulsiveChoosers_BisBas)),basD(~impulsiveChoosers_BisBas),10,rgb('orangered'),'filled')
    scatter(2*ones(1,sum(impulsiveChoosers_BisBas)),basF(impulsiveChoosers_BisBas),10,rgb('red'),'filled','^')
    scatter(2.3*ones(1,sum(~impulsiveChoosers_BisBas)),basF(~impulsiveChoosers_BisBas),10,rgb('red'),'filled')
    scatter(3*ones(1,sum(impulsiveChoosers_BisBas)),basR(impulsiveChoosers_BisBas),10,rgb('black'),'filled','^')
    scatter(3.3*ones(1,sum(~impulsiveChoosers_BisBas)),basR(~impulsiveChoosers_BisBas),10,rgb('black'),'filled')
    scatter(4*ones(1,sum(impulsiveChoosers_BisBas)),bisScore(impulsiveChoosers_BisBas),10,rgb('green'),'filled','^')
    scatter(4.3*ones(1,sum(~impulsiveChoosers_BisBas)),bisScore(~impulsiveChoosers_BisBas),10,rgb('green'),'filled')
end
hold off
xticks([1 2 3 4])
set(gca,'XTickLabel',{'BasDrive','BasFun','BasReward','BisScore'})
xlim([0.5 5])
ylabel('BISBAS SCORES')
ylim([0 30])
title('BISBAS subscale scores by impulsivity')
% BISBAS scores grouped by impulsivtiy metrics

% BISBAS total Score grouped by impulsivity metrics
subplot(2,3,[3 6])
hold on
if BOXPLOT
    betterBoxplot(1,totalScore(impulsiveChoosers_BisBas),rgb('black'),15,'^',1,false)
    betterBoxplot(1.5,totalScore(~impulsiveChoosers_BisBas),rgb('black'),15,'o',1,false)
else
    scatter(ones(1,sum(impulsiveChoosers_BisBas)),totalScore(impulsiveChoosers_BisBas),10,rgb('gold'),'filled','^')
    scatter(1.3*ones(1,sum(~impulsiveChoosers_BisBas)),totalScore(~impulsiveChoosers_BisBas),10,rgb('gold'),'filled')
end
hold off
xticks([1.25])
set(gca,'XTickLabel',{'Total Scores'})
xlim([.75 1.75])
ylabel('TOTAL SCORE')
ylim([0 100])
title('BISBAS total scores by impulsivity')

%saveas(1,'D:\Data\Rhiannon\BART_RLDM_outputs\BISBAS\boxplotsBisBasImpulsivity_allpts.pdf')

%% Regress impulsivity scores with BISBAS scores %%
% BAS DRIVE
fprintf('\ntotal Bas Drive score vs. impulsivty metric\n')
lm = fitlm(impulsivityMetric_BisBas,basD)
figure(2)
subplot(3,3,1)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('impulsivity metric')
ylabel('Bas Drive score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% BAS FUN SEEKING
fprintf('\ntotal Bas Fun score vs. impulsivty metric\n')
lm = fitlm(impulsivityMetric_BisBas,basF)
figure(2)
subplot(3,3,2)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('impulsivity metric')
ylabel('Bas Fun Seeking score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% BAS REWARD RESPONSIVENESS
fprintf('\ntotal Bas Reward score vs. impulsivty metric\n')
lm = fitlm(impulsivityMetric_BisBas,basR)
figure(2)
subplot(3,3,3)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('impulsivity metric')
ylabel('Bas Reward Responsiveness score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% BIS SCORES
fprintf('\ntotal Bis score vs. impulsivty metric\n')
lm = fitlm(impulsivityMetric_BisBas,bisScore)
figure(2)
subplot(3,3,4)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('impulsivity metric')
ylabel('Bis score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% TOTAL SCORES
fprintf('\ntotal BISBAS score vs. impulsivty metric\n')
lm = fitlm(impulsivityMetric_BisBas,totalScore)
figure(2)
subplot(3,3,5)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('impulsivity metric')
ylabel('total score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

%saveas(2,'D:\Data\Rhiannon\BART_RLDM_outputs\BISBAS\lmPlotsBisBasImpulsivity_allpts.pdf')

%% Regress BART metrics with BISBAS total scores %%

%Indexing for Plots
freePoints = [bisbasBhvStruct.freePoints];
freeITs = [bisbasBhvStruct.freeITs];
Accuracy = [bisbasBhvStruct.accuracyTot];
redAccuracy = [bisbasBhvStruct.redAccuracy];

% Free Points
fprintf('\nBISBAS total score vs. free points\n')
lm = fitlm(totalScore,freePoints)
figure(3)
subplot(3,3,1)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('free points')
ylabel('BISBAS total score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% free ITs
fprintf('\nBISBAS total score vs. free inflation times\n')
lm = fitlm(totalScore,freeITs)
figure(3)
subplot(3,3,2)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('free ITs')
ylabel('BISBAS total score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% Accuracy
fprintf('\nBISBAS total score vs. accuracy\n')
lm = fitlm(totalScore,Accuracy)
figure(3)
subplot(3,3,3)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('Accuracy')
ylabel('BISBAS total score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% Red Accuracy
fprintf('\nBISBAS total score vs. red accuracy\n')
lm = fitlm(totalScore,redAccuracy)
figure(3)
subplot(3,3,4)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('Red Accuracy')
ylabel('BISBAS total score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

%saveas(3,'D:\Data\Rhiannon\BART_RLDM_outputs\BISBAS\lmPlotsBisBasTotalScore_BARTbehavior_allpts.pdf')

%% Regress BART metrics with BISBAS total scores %%

% Free Points
fprintf('\nBISBAS Reward score vs. free points\n')
lm = fitlm(basR,freePoints)
figure(4)
subplot(3,3,1)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('free points')
ylabel('BISBAS Reward score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% free ITs
fprintf('\nBISBAS Reward score vs. free inflation times\n')
lm = fitlm(basR,freeITs)
figure(4)
subplot(3,3,2)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('free ITs')
ylabel('BISBAS Reward score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

% Accuracy
fprintf('\nBISBAS Reward score vs. accuracy\n')
lm = fitlm(basR,Accuracy)
figure(4)
subplot(3,3,3)
h = lm.plot;
h(1).Marker = '.';
h(1).Color = 'k';
xlabel('Accuracy')
ylabel('BISBAS Reward score')
legend off
title(sprintf('p = %.4f',lm.Coefficients{"x1","pValue"}))
axis tight square

keyboard

%saveas(4,'D:\Data\Rhiannon\BART_RLDM_outputs\BISBAS\lmPlotsBisBasRewardScore_BARTbehavior_allpts.pdf')


%% regressions for accuracy/points for BISBAS scores
%% regressions for color points for BISBAS scores