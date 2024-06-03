% function [] = compareRTs(nevFile,matFile)

% ::blackrock data "channel" codes::
%1: trial start ::      [1 2 3 4 11 12 13 14] = [Y O R G Yc Oc Rc Gc]
%2: responded ::        [22]
%3: inflating ::        [23 24] = [start stop]
%4: banked ::           [25]
%5: popped ::           [26]
%6: outcome shown ::    [100 101] = [correct incorrect] 
%7: max rt exceeded  :: [127]
%8: trial over ::		[120]

% loading nevFile
ptID = '201810';
nevFile = '/media/user1/data4TB/data/BART/BART_EMU/patient1/datafile_real001.nev';
NEV = openNEV(nevFile);

% loading matFile
matFile = '/media/user1/data4TB/data/BART/BART_EMU/patient1/test.2.bartc.mat';
load(matFile)

% digital triggers.
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% [20170713] I made a small error in the handle_input2.m script such that
% sometimes there is an additional 23 after the initial 23 (signifying the
% start of the inflation). The next two lines remove that second 23. 
infIdx = trigs==23;
infIdx([diff(infIdx)==0; false] & infIdx==1) = 0;

% 
respTimes = trigTimes(trigs==25 | trigs==26);
trialStarts = trigTimes(trigs<15);
trialStarts = trialStarts(2:length(respTimes)+1);
inflateStarts = trigTimes(infIdx);


% calculating RTs. 
brRTs = respTimes-trialStarts;
mlRTs = [data(1:end-1).inflate_time];
ctrlIdcs = ~logical([data(1:end-1).is_control]);

maxRTerror = max(abs(brRTs(ctrlIdcs)-mlRTs(ctrlIdcs)));

% plotting RTs. 
figure
hold on
scatter(brRTs(ctrlIdcs),mlRTs(ctrlIdcs))
plot(1:7,1:7,'--k')
hold off
xlabel('RTs calculated from blackrock')
ylabel('RTs in mat data')



