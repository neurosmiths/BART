% validate task structure [20181212]
%% start by loading task data
ptID = '201811';
nevList = dir(sprintf('~/data/BART/BART_EMU/%s/Data/*.nev',ptID));
if length(nevList)>1
    error('many nev files available for this patient. Please specify...')
elseif length(nevList)<1
    error('no nev files found...')
else
    nevFile = fullfile(nevList.folder,nevList.name);
end


%% parsing triggers
NEV = openNEV(nevFile,'overwrite');
if NEV.Data.SerialDigitalIO.UnparsedData(end)<15
    trigs = NEV.Data.SerialDigitalIO.UnparsedData(1:end-1);
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec(1:end-1);
else
    trigs = NEV.Data.SerialDigitalIO.UnparsedData;
    trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;
end


if ishandle(1); close(1); end
%% plotting triggers for the whole session.
trigTypes = unique(trigs);
for tr = 1:length(trigTypes)
    switch trigTypes(tr)
        case 0
            trCol = rgb('white');
            trName = '';
        case 1
            trCol = rgb('gold');
            trName = 'yellow';
        case 2
            trCol = rgb('orangered');
            trName = 'orange';
        case 3
            trCol = rgb('red');
            trName = 'red';
        case 11
            trCol = rgb('khaki');
            trName = 'yellowCTRL';
        case 12
            trCol = rgb('coral');
            trName = 'orangeCTRL';
        case 13
            trCol = rgb('lightcoral');
            trName = 'redCTRL';
        case 14
            trCol = rgb('dimgray');
            trName = 'grayCTRL';
        case 23
            trCol = rgb('deepPink');
            trName = 'inflateSTART (23)';
        case 24
            trCol = rgb('darkMagenta');
            trName = 'inflateSTOP (24)';
        case 25
            trCol = rgb('lime');
            trName = 'banked (25)';
        case 26
            trCol = rgb('tomato');
            trName = 'popped (26)';
        case 100
            trCol = rgb('darkseagreen');
            trName = 'correctOutcome (100)';
        case 101
            trCol = rgb('rosybrown');
            trName = 'incorrectOutcome (101)';
        case 120
            trCol = rgb('black');
            trName = 'trialEND (120)';
    end
    %% current trigger
    thisTrig = trigTypes(tr);
    
    
    figure(1)
    %% plotting all triggers for the full session
    subplot(2,1,1)
    hold on
    trTimes = trigTimes(trigs==thisTrig);
    for lb = 1:length(trTimes)
        line([trTimes(lb) trTimes(lb)],[0 1],'color',trCol,'linewidth',1)
        h=text(trTimes(lb),1,trName,'color',trCol,'fontsize',8);
        set(h,'Rotation',45);
    end
    hold off
    axis tight
    xlabel('session timeline (s)')
    ylim([0 2])
    ylabel('marker IDs')
    
end


%% setting up variables relative to the appearance of balloons.
% [20170713] I made a small error in the handle_input2.m script such that
% sometimes there is an additional 23 after the initial 23 (signifying the
% start of the inflation). The next two lines remove that second 23.
infIdx = trigs==23;
infIdx([diff(infIdx)==0; false] & infIdx==1) = 0;

% response times
if strcmp(ptID(1:4),'2018')
    respTimes = trigTimes(trigs==25 | trigs==26);
else
%     respTimes = 
end    
trialStarts = trigTimes(trigs<15);
trialStarts = trialStarts(2:length(respTimes)+1);
inflateStarts = trigTimes(infIdx);
nTrials = length(respTimes);

% calculating RTs.
brRTs = respTimes-trialStarts;

% [20181213] it appears that the inflate stop markers atually appear at the
% time the baloons stop inflating.
yellowITs = trigTimes(find(trigs==1)+1) - trigTimes(trigs==1);
orangeITs = trigTimes(find(trigs==2)+1) - trigTimes(trigs==2);
redITs = trigTimes(find(trigs==3)+1) - trigTimes(trigs==3);
yellowCTRLITs = trigTimes(find(trigs==11)+1) - trigTimes(trigs==11);
orangeCTRLITs = trigTimes(find(trigs==12)+1) - trigTimes(trigs==12);
redCTRLITs = trigTimes(find(trigs==13)+1) - trigTimes(trigs==13);
grayCTRLITs = trigTimes(find(trigs==14)+1) - trigTimes(trigs==14);

% histogram params
nBins = 50;
histoAlpha = 0.5;

% plotting mean timing for each trial
subplot(2,1,2)
hold on
yith = histfit(yellowITs,nBins);
oith = histfit(orangeITs,nBins);
rith = histfit(redITs,nBins);
ycith = histfit(yellowCTRLITs,nBins);
ocith = histfit(orangeCTRLITs,nBins);
rcith = histfit(redCTRLITs,nBins);
gcith = histfit(grayCTRLITs,nBins);
hold off

% histo deets
rith(1).FaceColor = rgb('red');
rith(1).FaceAlpha = histoAlpha;
rith(1).EdgeColor = 'none';
oith(1).FaceColor = rgb('orangered');
oith(1).FaceAlpha = histoAlpha;
oith(1).EdgeColor = 'none';
yith(1).FaceColor = rgb('gold');
yith(1).FaceAlpha = histoAlpha;
yith(1).EdgeColor = 'none';
% line deets
rith(2).Color = rgb('red');
rith(2).LineWidth = 1;
oith(2).Color = rgb('orangered');
oith(2).LineWidth = 1;
yith(2).Color = rgb('gold');

% control histo deets
rcith(1).FaceColor = rgb('lightcoral');
rcith(1).FaceAlpha = histoAlpha;
rcith(1).EdgeColor = 'none';
ocith(1).FaceColor = rgb('coral');
ocith(1).FaceAlpha = histoAlpha;
ocith(1).EdgeColor = 'none';
ycith(1).FaceColor = rgb('khaki');
ycith(1).FaceAlpha = histoAlpha;
ycith(1).EdgeColor = 'none';
% control line deets
rcith(2).Color = rgb('lightcoral');
rcith(2).LineWidth = 1;
ocith(2).Color = rgb('coral');
ocith(2).LineWidth = 1;
ycith(2).Color = rgb('khaki');
% gray deets
gcith(1).FaceColor = rgb('dimgray');
gcith(1).FaceAlpha = histoAlpha;
gcith(1).EdgeColor = 'none';
gcith(2).Color = rgb('dimgray');
gcith(2).LineWidth = 1;

% axis deets
xlabel('inflate time (s)')
ylabel('count')

maximize(1)