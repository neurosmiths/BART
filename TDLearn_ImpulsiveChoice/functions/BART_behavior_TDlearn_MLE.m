function [TDdata] = BART_behavior_TDlearn_MLE(ptID,whichTD)
% BART_BEHAVIOR_TDLEARN analyzes BART behavior from the perspective of a 
%	temporal difference learning RL algorithm. The second input argument
%	is a string indicating whetehr to do 'vanilla' TD learning or 
%	'distributional' TD learning. 

% author: EHS20201023


%% extra info about task beahvioral markers. 
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

% patient details.
% ptID = '201915';

nevList = dir(sprintf('~/data/BART/BART_EMU/%s/Data/*.nev',ptID));
nevFile = fullfile(nevList.folder,nevList.name);
% trodeLabels = ptTrodesBART(ptID);

% initializing bhv output
TDdata.patientID = ptID;
TDdata.type = whichTD;

% loading matFile
matFile = sprintf('~/data/BART/BART_EMU/%s/Data/%s.bartBHV.mat',ptID,ptID);
load(matFile)

% load and define triggers from nevFle
NEV = openNEV(nevFile,'overwrite');
trigs = NEV.Data.SerialDigitalIO.UnparsedData;
trigTimes = NEV.Data.SerialDigitalIO.TimeStampSec;

% [20170713] I made a small error in the handle_input2.m script such that
% sometimes there is an additional 23 after the initial 23 (signifying the
% start of the inflation). The next two lines remove that second 23.
infIdx = trigs==23;
infIdx([diff(infIdx)==0; false] & infIdx==1) = 0;

% trial IDs
isCTRL = logical([data.is_control]);
balloonIDs = trigs(trigs==1 | trigs==2 | trigs==3 | trigs==11 | trigs==12 | trigs==13 | trigs==14);
outcomeType = trigs(sort([find(trigs==25); find(trigs==26)]))-24;

% response times
respTimes = trigTimes(trigs==25 | trigs==26);
trialStarts = trigTimes(trigs<15);
trialStarts = trialStarts(2:length(respTimes)+1);
inflateStarts = trigTimes(infIdx);
nTrials = length(respTimes);

% adjusting for numbers of trials. 
balloonIDs = balloonIDs(1:nTrials);
isCTRL = isCTRL(1:nTrials);

% risk colormap
cMap(2,:) = [1 0.9 0];
cMap(3,:) = [1 0.5 0];
cMap(4,:) = [1 0 0];
cMap(1,:) = [0.5 0.5 0.5];

% risk and reward colormap
rcMap(1,:) = [0.5 0.5 0.5];
rcMap(2,:) = [1 0 0];
rcMap(3,:) = [1 0.5 0];
rcMap(4,:) = [1 0.9 0];
rcMap(5,:) = rgb('lightcoral');
rcMap(6,:) = rgb('rosybrown');
rcMap(7,:) = rgb('violet');

% colormap per trial. 
balloonEdgeColorMap = ones(length(balloonIDs),3)*.5;
riskRewardEdgeColorMap = zeros(length(balloonIDs),3);

% populating balloon color map
for x = 1:3
	balloonEdgeColorMap(balloonIDs==x,:) = repmat(cMap(x+1,:),sum(balloonIDs==x),1);
end

%% trying to measure optimized learning rates with MLE
R = data.points;


RewardPE(t) = R(t) - V(t-1);





% fitting TD learning models.
if strcmp(whichTD,'vanilla')
		% Initalize variables. 
		R(1) = data(1).points;
		X(1) = 0.5;								% starting with a coin flip for risk 
		V = zeros(1,nTrials);
		Vrisk = zeros(1,nTrials);
		a = 0.1;

		for t = 2:nTrials
			% reward varriable. 
			if strcmp(data(t).result,'banked')
				R(t) = data(t).points;			% outcome on current trial. 
			else
				R(t) = 0;
			end

			% risk variable. 
			if isCTRL(t)
				X(t) = 0;
			else
				X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
				% current trial risk is defined as P(pop) on previously observed balloons. 
			end
			
			% updating risk and reward PE
			RewardPE(t) = R(t) - V(t-1);
			XPE(t) = X(t) - Vrisk(t-1);

			% updating valua
			Vrisk(t) = Vrisk(t-1) + a*XPE(t-1);
			V(t) = V(t-1) + a*RewardPE(t-1);


			% if t==nTrials; keyboard; end;
		end

		% stats
		PElm = fitlm(RewardPE(50:end)',XPE(50:end)')	% This data appears to cluster rather than fit a line, so going to plot with balloon color
		% [20201029]:: if the clusters look neat, it may be worth looking at separability vs. impulsivity. 
		Vlm = fitlm(V(50:end)',Vrisk(50:end)')

		% plotting 
		figure
		subplot(3,6,[1 2])
		plot(1:nTrials,RewardPE)
		yl = ylim;
		xlabel('trials')
		ylabel('reward PE')
		title(sprintf('patient: %s',ptID))

		subplot(3,6,3)
		histogram(RewardPE,20,'Orientation','horizontal','DisplayStyle','stairs')
		ylim(yl);

		subplot(3,6,[4 5])
		plot(1:nTrials,V)
		yl = ylim;
		xlabel('trials')
		ylabel('value estimate')

		subplot(3,6,6)
		histogram(V,20,'Orientation','horizontal','DisplayStyle','stairs')
		ylim(yl);

		subplot(3,6,[7 8])
		plot(1:nTrials,XPE)
		yl = ylim;
		xlabel('trials')
		ylabel('Risk PE')

		subplot(3,6,9)
		histogram(XPE,20,'Orientation','horizontal','DisplayStyle','stairs')
		ylim(yl);

		subplot(3,6,[10 11])
		plot(1:nTrials,Vrisk)
		yl = ylim;
		xlabel('trials')
		ylabel('risk estimate')

		subplot(3,6,12)
		histogram(Vrisk,20,'Orientation','horizontal','DisplayStyle','stairs')
		ylim(yl);

		subplot(3,2,5)
		hold on
		scatter(RewardPE(outcomeType==1),XPE(outcomeType==1),8,balloonEdgeColorMap(outcomeType==1,:),'filled')
		scatter(RewardPE(outcomeType==2),XPE(outcomeType==2),10,balloonEdgeColorMap(outcomeType==2,:))
		hold off
		xlabel('reward PE')
		ylabel('risk PE')
		title('risk and reward PEs, colored by trial balloon color')
		%title(sprintf('t(%d) = %.2f, p = %.2f, R^2 = %.2f',PElm.DFE,PElm.Coefficients{2,3},PElm.Coefficients{2,4},PElm.Rsquared.Ordinary))

		subplot(3,2,6)
		Vlm.plot
		legend('off')
		xlabel('Vreward')
		ylabel('Vrisk')
		title(sprintf('t(%d) = %.2f, p = %.2f, R^2 = %.2f',Vlm.DFE,Vlm.Coefficients{2,3},Vlm.Coefficients{2,4},Vlm.Rsquared.Ordinary))
	
		% saving figures
		halfMaximize(gcf,'left')
		saveas(gcf,sprintf('~/Dropbox/BART_distRLbehavior/vanilla//%s_simpleTDRL.pdf',ptID))

		% saving variables in a struct for output
		TDdata.nTrials = nTrials;
		TDdata.a = a;
		TDdata.R = R;
		TDdata.V = V;
		TDdata.RewardPE = RewardPE;
		TDdata.perTrialRisk = X;
		TDdata.RiskPE = XPE;
		TDdata.Vrisk = Vrisk; 
		TDdata.rewardVsRiskPElm = PElm;
		TDdata.valueVsRiskEstimateLM = Vlm;


elseif strcmp(whichTD,'distributional')
		% expectile or quantile learning rule
		whichDistRL = 'expectile';
	
		% Initalize variables. 
		R(1) = data(1).points
		X(1) = 0.5
		V = zeros(1,nTrials);
		Vrisk = zeros(1,nTrials);
		
		% looping over trials. 
		for t = 2:nTrials
			if strcmp(data(t).result,'banked')
				R(t) = data(t).points;	% outcome on current trial. 
			else
				R(t) = 0;						%TODO:: replace this with negative value (points lost)
			end
	
			% risk variable. 
			if isCTRL(t)
				X(t) = 0;
			else
				X(t) = sum(outcomeType(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10)==1)./sum(balloonIDs==balloonIDs(t) | balloonIDs==balloonIDs(t)+10);
				% current trial risk is defined as P(pop) on previously observed balloons. 
			end

			% update PEs
			RewardPE(t) = R(t) - V(t-1);
			XPE(t) = X(t) - Vrisk(t-1);

			% specifying overall alphas. 
			alphas = [0.1 0.9]';

			% updating values for risk. 
			if strcmp(whichDistRL,'quantile')
				if XPE(t)<=0
					Vrisk(t) = Vrisk(t-1) - alphas(1);
				elseif XPE(t)>0
					Vrisk(t) = Vrisk(t-1) + alphas(2);
				end
			elseif strcmp(whichDistRL,'expectile')
				if XPE(t)<=0
						Vrisk(t) = Vrisk(t-1) + alphas(1).*XPE(t-1);
				elseif XPE(t)>0
						Vrisk(t) = Vrisk(t-1) + alphas(2).*XPE(t-1);
				end
			end
			
			% updating values for reward. 
			if strcmp(whichDistRL,'quantile')
				if RewardPE(t)<=0
					V(t) = V(t-1) - alphas(1);
				elseif RewardPE(t)>0
					V(t) = V(t-1) + alphas(2);
				end
			elseif strcmp(whichDistRL,'expectile')
				if RewardPE(t)<=0
						V(t) = V(t-1) + alphas(1).*RewardPE(t-1);
				elseif RewardPE(t)>0
						V(t) = V(t-1) + alphas(2).*RewardPE(t-1);
				end
			end
		end

		% plotting
		figure
		subplot(2,3,1)
		plot(1:nTrials,RewardPE);
		xlabel('trials')
		ylabel('reward PE')
		title(sprintf('patient: %s',ptID));

		subplot(2,3,2)
		hold on
		plot(1:nTrials,V);
		xlabel('trials')
		ylabel('value estimate')
		title(sprintf('%s learning rule, alphas: %.2f',whichDistRL,alphas))

		subplot(2,3,4)
		plot(1:nTrials,XPE);
		xlabel('trials')
		ylabel('Risk PE')

		subplot(2,3,5)
		plot(1:nTrials,Vrisk);
		xlabel('trials')
		ylabel('risk estimate')

 		subplot(2,6,5)
		hold on
		h = histogram(RewardPE,'DisplayStyle','stairs')
		h.Orientation = 'horizontal';
		hold off
		title('RPE')

 		subplot(2,6,6)
		hold on
		h = histogram(V,'DisplayStyle','stairs')
		h.Orientation = 'horizontal';
		hold off
		title('Value')

 		subplot(2,6,11)
		hold on
		h = histogram(XPE,'DisplayStyle','stairs')
		h.Orientation = 'horizontal';
		hold off
		title('risk PE')

 		subplot(2,6,12)
		hold on
		h = histogram(Vrisk,'DisplayStyle','stairs')
		h.Orientation = 'horizontal';
		hold off
		title('Risk')

		% saving figures
		orient(gcf,'landscape')
		saveas(gcf,sprintf('~/Dropbox/BART_distRLbehavior/%s_distributionalTDRL_%s.pdf',ptID,whichDistRL))

		% saving variables in a struct for output
		TDdata.nTrials = nTrials;
		TDdata.a = alphas
		TDdata.R = R;
		TDdata.V = V;
		TDdata.RewardPE = RewardPE;
		TDdata.perTrialRisk = X;
		TDdata.RiskPE = XPE;
		TDdata.Vrisk = Vrisk; 
	%	TDdata.rewardVsRiskPElm = PElm;
	%	TDdata.valueVsRiskEstimateLM = Vlm;

end




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%~~ OLD CODE ~~%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % calculating RTs.
% brRTs = respTimes-trialStarts;
% mlRTs = [data(1:nTrials).inflate_time];
% ctrlIdcs = ~logical([data(1:nTrials).is_control]);
% 
% maxRTerror = max(abs(brRTs(ctrlIdcs)-mlRTs(ctrlIdcs)));
% 
% 
% % plotting differences in RTs determined from the blackrock and mat data
% compareRTs = false;
% if compareRTs
%     figure(678)
%     hold on
%     scatter(brRTs(ctrlIdcs),mlRTs(ctrlIdcs))
%     plot(1:7,1:7,'--k')
%     hold off
%     xlabel('RTs calculated from blackrock')
%     ylabel('RTs in mat data')
%     pause(30)
%     close(678)
% end
% 
% 
% 
% %% [20181009] first organize data
% % overall parameters
% figID = 1123;
% histoAlpha = 0.4;
% % in case we want to specify the number of histogram bins. 
% % uses these many bins for K-L Divergence impulsivity metric.
% % leave empty to use the defaults from matlab. ~10-12 is reasonable.
% nBins = []; 
% 
% % caluclating results.
% data = data(1:nTrials);
% isCtrl = logical([data.is_control]);
% nCtrls = sum(isCtrl);
% trialTypeWithCtrls = [data.trial_type];
% % 3 = red; 2 = orange; 1 = yellow; 4 = control
% 
% % 1) accuracy
% banked = double(strcmp({data.result},'banked'));
% popped = double(strcmp({data.result},'popped'));
% nPopped = sum(popped)
% 
% % overall accuracy percentages
% redAccuracy = (sum(banked(trialTypeWithCtrls==3))./sum(trialTypeWithCtrls==3))*100;
% orangeAccuracy = (sum(banked(trialTypeWithCtrls==2))./sum(trialTypeWithCtrls==2))*100;
% yellowAccuracy =  (sum(banked(trialTypeWithCtrls==1))./sum(trialTypeWithCtrls==1))*100;
% 
% % 2) RT
% allRTs = [data.rt];
% freeRTs = allRTs(~isCtrl);
% ctrlRTs = allRTs(isCtrl);
% redRTs = allRTs(trialTypeWithCtrls==3);
% orangeRTs = allRTs(trialTypeWithCtrls==2);
% yellowRTs = allRTs(trialTypeWithCtrls==1);
% % [20181009] RT isn't really that interesting across conditions. maybe
% % just report as text...
% if logical(banked(end))
%     postBankRTs = allRTs(find(banked(1:end-1))+1);
%     postPopRTs = allRTs(find(popped)+1);
% else
%     postBankRTs = allRTs(find(banked)+1);
%     postPopRTs = allRTs(find(popped(1:end-1))+1);
% end
% 
% % 3) IT
% allITs = [data.inflate_time];
% freeITs = allITs(~isCtrl);
% ctrlITs = allITs(isCtrl);
% redITs = allITs(trialTypeWithCtrls==3 & ~isCtrl);
% orangeITs = allITs(trialTypeWithCtrls==2 & ~isCtrl);
% yellowITs = allITs(trialTypeWithCtrls==1 & ~isCtrl);
% if logical(banked(end))
%     postBankITs = allITs(find(banked(1:end-1))+1);
%     postPopITs = allITs(find(popped)+1);
% else
%     postBankITs = allITs(find(banked)+1);
%     postPopITs = allITs(find(popped(1:end-1))+1);
% end
% 
% % 4) points (proxy for IT)
% allPoints = [data.points];
% bhvStruct.allPoints = sum(allPoints);
% bhvStruct.freePoints = sum(allPoints(~isCtrl));
% bhvStruct.ctrlPoints = sum(allPoints(isCtrl));
% bhvStruct.redPoints = allPoints(trialTypeWithCtrls==3);
% bhvStruct.orangePoints = allPoints(trialTypeWithCtrls==2);
% bhvStruct.yellowPoints = allPoints(trialTypeWithCtrls==1);
% 
% 
% %% [20190311] adding values to the bhvStruct - more later...
% bhvStruct.totalTrials = nTrials;
% bhvStruct.nCtrls = nCtrls;
% bhvStruct.nActiveTrials = nTrials/nCtrls;
% bhvStruct.redAccuracy = redAccuracy;
% bhvStruct.orangeAccuracy = orangeAccuracy;
% bhvStruct.yellowAccuracy = yellowAccuracy;
% 
% 
% %% more stats
% % histograms for kldiv
% maxIT = max(allITs);
% nBns = 10;
% [ctrlHist,edges] = histcounts(ctrlITs,0:maxIT./nBns:maxIT);
% [freeHist,edges] = histcounts(freeITs,0:maxIT./nBns:maxIT);
% 
% c = midpoint(edges);
% centers = edges(1:end-1)+c(1);
% eps = 0.00000000000000000001;
% d=sum((freeHist./sum(freeHist)).*log2(freeHist./sum(freeHist)+eps)-(freeHist./sum(freeHist)).*log2(ctrlHist./sum(ctrlHist)+eps)); %KL(h1,h2)
% bhvStruct.impulsivityKLD = d; % *(sum(banked)./nTrials)
% [ITp,~,ITstats] = ranksum(ctrlITs,freeITs);
% bhvStruct.activeVsPassiveStats.pVal = ITp;
% bhvStruct.activeVsPassiveStats.zVal = ITstats.zval;
% 
% 
% %% plot a summary figure. 
% % if the figure exists, close it. 
% if ishandle(figID); close(figID); end
% figure(figID)
% 
% %% text in the top right square
% plotmultipleaxes(1,2,3,0.05,figID)
% % stats
% redBanked = zeros(1,length(banked));
% redBanked(trialTypeWithCtrls==3)=1;
% redBanked(logical(popped))=0;
% orangeBanked = zeros(1,length(banked));
% orangeBanked(trialTypeWithCtrls==2)=1;
% orangeBanked(logical(popped))=0;
% yellowBanked = zeros(1,length(banked));
% yellowBanked(trialTypeWithCtrls==1)=1;
% yellowBanked(logical(popped))=0;
% [accuracyContingencyTbl,chi2Accuracy,pAccuracy] = crosstab(redBanked,orangeBanked,yellowBanked);
% % plotted text
% title('BART behavior:')
% hold on
% maxLim = 24;
% text(0,maxLim,sprintf('Subject: %s',ptID));
% text(0,maxLim-2,sprintf('total baloons = %d (%d controls)',nTrials,nCtrls));
% text(0,maxLim-4,sprintf('Accuracy for red balloons: %.1f percent.',redAccuracy));
% text(0,maxLim-6,sprintf('Accuracy for orange balloons: %.1f percent.',orangeAccuracy));
% text(0,maxLim-8,sprintf('Accuracy for yellow balloons: %.1f percent.',yellowAccuracy));
% if pAccuracy<0.01
%     text(0,maxLim-10,sprintf('Signifcant difference in accuracy among balloon colors: X^2(%d)=%.1f, p=%.2d',nTrials-2,chi2Accuracy,pAccuracy));
% else
%     text(0,maxLim-10,sprintf('NO signifcant difference among balloon colors: X^2(%d)=%.1f, p=%.2d',nTrials-2,chi2Accuracy,pAccuracy));
% end
% text(0,maxLim-12,sprintf('Total accuracy: %.1f percent.',(sum(banked)./nTrials)*100));
% text(0,maxLim-14,sprintf('KLD-derived impulsivity = %.2f ',d));
% hold off
% % deets
% ylim([0 25])
% axis off
% 
% %% adding stats to the bhvStruct
% bhvStruct.accuracyTot = (sum(banked)./nTrials)*100;
% bhvStruct.accuracyStats.df = nTrials-2;
% bhvStruct.accuracyStats.chi2Acciracy = chi2Accuracy;
% bhvStruct.accuracyStats.pAccuracy = pAccuracy;
% 
% 
% %% inflate time histogramsfor free vs. ctrls
% plotmultipleaxes(2,2,3,0.05,figID)
% % stats
% title(sprintf('ranksum test, U=%d, p=%.2d.',ITstats.ranksum,ITp))
% % plots
% hold on
% cith = histfit(ctrlITs,nBins);
% fith = histfit(freeITs,nBins);
% 
% hold off
% % deets
% cith(1).FaceEdgeColor = rgb('forestgreen');
% cith(1).FaceAlpha = histoAlpha;
% cith(2).EdgeColor = rgb('forestgreen');
% cith(2).LineWidth = 1;
% fith(1).FaceEdgeColor = rgb('navy');
% fith(1).FaceAlpha = histoAlpha;
% fith(2).EdgeColor = rgb('navy');
% fith(2).LineWidth = 1;
% axis square
% ylabel('IT count')
% xlabel('time (s)')
% xlim([0 max(allITs)])
% 
% 
% %% inflate time histograms for baloon types
% plotmultipleaxes(5,2,3,0.05,figID)
% % stats
% [pColIT,tblColIT,statsColIT] = kruskalwallis(allITs(trialTypeWithCtrls<4 & ~isCtrl),trialTypeWithCtrls(trialTypeWithCtrls<4 & ~isCtrl),'off');
% title(sprintf('excluding passive, kruskal-wallis, X^2=%.2f, p=%.2d.',tblColIT{2,5},tblColIT{2,6}))
% % plots
% hold on
% yith = histfit(yellowITs,nBins);
% oith = histfit(orangeITs,nBins);
% rith = histfit(redITs,nBins);
% hold off
% % deets
% rith(1).FaceEdgeColor = rgb('red');
% rith(1).FaceAlpha = histoAlpha;
% oith(1).FaceEdgeColor = rgb('orangered');
% oith(1).FaceAlpha = histoAlpha;
% yith(1).FaceEdgeColor = rgb('gold');
% yith(1).FaceAlpha = histoAlpha;
% rith(2).EdgeColor = rgb('red');
% rith(2).LineWidth = 1;
% oith(2).EdgeColor = rgb('orangered');
% oith(2).LineWidth = 1;
% yith(2).EdgeColor = rgb('gold');
% yith(2).LineWidth = 1;
% axis square
% ylabel('IT count')
% xlabel('time (s)')
% xlim([0 max(allITs)])
% 
% 
% %% more stats
% bhvStruct.balloonEdgeColorInflationTimeStats.pVal = pColIT;
% bhvStruct.balloonEdgeColorInflationTimeStats.statTbl = tblColIT;
% bhvStruct.balloonEdgeColorInflationTimeStats.statStruct = statsColIT;
% 
% 
% %% [20181009] plot [RT] distributions as violin plots
% plotmultipleaxes(15,7,3,0.05,figID)
% % stats
% [RTp,~,RTstats] = ranksum(ctrlRTs,freeRTs);
% % plots
% hold on
% violinPlot(ctrlRTs','showMM',4,'color',rgb('forestgreen'),'xValues',1)
% violinPlot(freeRTs','showMM',4,'color',rgb('navy'),'xValues',2)
% hold off
% % deets
% ylabel('RT (s)')
% xlabel(sprintf('ranksum, U(2)=%d, p=%.2f.',RTstats.ranksum,RTp))
% if ceil(max(allRTs))>10
%     ylim([0 10])
% else
%     ylim([0 ceil(max(allRTs))])
% end
% 
% 
% %% more stats
% bhvStruct.reactionTimeStats.pVal = RTp;
% bhvStruct.reactionTimeSTats.zVal = RTstats.zval;
% 
% 
% %% [20190116] plot [RT] distributions as violin plots
% plotmultipleaxes(18,7,3,0.05,figID)
% % stats
% [pColRT,tblColRT,statsColRT] = kruskalwallis(allRTs(trialTypeWithCtrls<4 & ~isCtrl),trialTypeWithCtrls(trialTypeWithCtrls<4 & ~isCtrl),'off');
% % plots
% hold on
% violinPlot(redRTs','showMM',4,'color',rgb('red'),'xValues',1)
% violinPlot(orangeRTs','showMM',4,'color',rgb('orangered'),'xValues',2)
% violinPlot(yellowRTs','showMM',4,'color',rgb('gold'),'xValues',3)
% hold off
% % deets
% ylabel('RT (s)')
% xlabel(sprintf('kruskal-wallis, X^2(2)=%.2f, p=%.2f.',tblColRT{2,5},tblColRT{2,6}))
% if ceil(max(allRTs))>10
%     ylim([0 10])
% else
%     ylim([0 ceil(max(allRTs))])
% end
% 
% 
% %% more stats
% bhvStruct.balloonEdgeColorRTStats.pVal = pColRT;
% bhvStruct.balloonEdgeColorRTStats.statStruct = statsColRT;
% 
% 
% %% [20181009] plot follownig trial [RT] distributions as violin plots
% plotmultipleaxes(21,7,3,0.05,figID)
% % stats
% [RTp,~,RTstats] = ranksum(postBankRTs,postPopRTs);
% % plots
% hold on
% violinPlot(postBankRTs','showMM',4,'color',rgb('forestgreen'),'xValues',1)
% violinPlot(postPopRTs','showMM',4,'color',rgb('orangered'),'xValues',2)
% hold off
% % deets
% ylabel('RT (s)')
% xlabel(sprintf('ranksum, U(2)=%d, p=%.2f.',RTstats.ranksum,RTp))
% if ceil(max(allRTs))>10
%     ylim([0 10])
% else
%     ylim([0 ceil(max(allRTs))])
% end
% 
% 
% %% more stats
% bhvStruct.postBankPopRTStats.pVal = RTp;
% bhvStruct.postBankPopRTStats.zVal = RTstats.zval;
% bhvStruct.postBankPopRTStats.bankPopMeans = [mean(postBankRTs) mean(postPopRTs)];
% bhvStruct.postBankPopRTStats.bankPopSTDs = [std(postBankRTs) std(postPopRTs)];
% 
% 
% % plotting points histograms for baloon types
% plotmultipleaxes(3,2,3,0.05,figID)
% % stats
% [pColPts,tblColPts,statsColPTs] = kruskalwallis(allPoints(trialTypeWithCtrls<4),trialTypeWithCtrls(trialTypeWithCtrls<4),'off');
% title(sprintf('including passive, kruskal-wallis, X^2=%.2f, p=%.2d.',tblColPts{2,5},tblColPts{2,6}))
% % plots
% hold on
% yith = histfit(bhvStruct.yellowPoints,nBins);
% oith = histfit(bhvStruct.orangePoints,nBins);
% rith = histfit(bhvStruct.redPoints,nBins);
% hold off
% % deets
% rith(1).FaceEdgeColor = rgb('red');
% rith(1).FaceAlpha = histoAlpha;
% oith(1).FaceEdgeColor = rgb('orangered');
% oith(1).FaceAlpha = histoAlpha;
% yith(1).FaceEdgeColor = rgb('gold');
% yith(1).FaceAlpha = histoAlpha;
% rith(2).EdgeColor = rgb('red');
% rith(2).LineWidth = 1;
% oith(2).EdgeColor = rgb('orangered');
% oith(2).LineWidth = 1;
% yith(2).EdgeColor = rgb('gold');
% yith(2).LineWidth = 1;
% axis square
% ylabel('point bin count')
% xlabel('points')
% xlim([0 max(allPoints)])
% 
% 
% %% saving figures
% halfMaximize(figID,'left')
% saveDir = sprintf('/media/user1/data4TB/Figs/BART',ptID);
% fName = sprintf('%s_BARTbehavior',ptID);
% print(fullfile(saveDir,fName),'-dpdf','-fillpage')
% 
% 
% %% More behavioral modeling for each patient. 
% % - model current inflation time as a fucntion of 1) previous trial
% % inflation time,
% % build table
% 
% % trialTypeCategorical = categorical(trialTypeWithCtrls);
% % trialTypeCategorical(trialTypeWithCtrls==1)= 'A';
% % trialTypeCategorical(trialTypeWithCtrls==2)= 'B';
% % trialTypeCategorical(trialTypeWithCtrls==3)= 'C';
% % trialTypeCategorical(trialTypeWithCtrls==4)= 'D';
% 
% pointsTable_wCTRL = table([data.points]',[data(2:end).points NaN]',trialTypeWithCtrls',...
%     'VariableNames',{'points','previousPoints','balloonType'});
% pointsTable = pointsTable_wCTRL(~isCtrl,:);
% 
% bhvStruct.glmPP = fitglme(pointsTable,'points ~ 1 + previousPoints + (1 + previousPoints|balloonType)') %  + (1 + 1|balloonType)
% bhvStruct.glmPP_p = bhvStruct.glmPP.Coefficients{2,6};
% bhvStruct.glmPP_t = bhvStruct.glmPP.Coefficients{2,4};
% bhvStruct.glmPP_R2 =  bhvStruct.glmPP.Rsquared.ordinary;
% 
% bhvStruct.glmBT = fitglme(pointsTable,'points ~ 1 + balloonType + (1 + 1|balloonType)')
% bhvStruct.glmBT_p = bhvStruct.glmBT.Coefficients{2,6};
% bhvStruct.glmBT_t = bhvStruct.glmBT.Coefficients{2,4};
% bhvStruct.glmBT_R2 =  bhvStruct.glmBT.Rsquared.ordinary;
% 
% bhvStruct.LMMcomparison = compare(bhvStruct.glmPP,bhvStruct.glmBT)
% 
% %% saving stats
% save(sprintf('/media/user1/data4TB/data/BART/BART_EMU/%s/Data/%s.bhvStruct.mat',ptID,ptID),'bhvStruct')
% 
% 
% 















