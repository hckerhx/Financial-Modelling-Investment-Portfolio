% Student Name: Lingkai Shen, Hangzuo Xiang, Haoran Gong
% Student ID: 1002066952, 1002376476 ,1002082836

clc
clear global
close all
format long
% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project2_Data_adjClose.csv', 'ReadVariableNames', true);
adjClose = table2timetable(adjClose);

factorRet = readtable('Project2_Data_FF_factors.csv', 'ReadVariableNames', true);
factorRet = table2timetable(factorRet);

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowTimes);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returnsT = array2table(returns);
returnsT.Properties.VariableNames = tickers;
returnsT = table2timetable(returnsT,'RowTimes',dates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of shares outstanding per company
mktNoShares = [36.40 9.86 15.05 10.91 45.44 15.43 43.45 33.26 51.57 45.25 ...
    69.38 8.58 64.56 30.80 36.92 7.70 3.32 54.68 38.99 5.78]' * 1e8;

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

methods={'MVO','robust','resampling','mostDiverse','CVaR'};
NoMethodsOld=size(methods,2);

% number of assets
n = size(tickers,1);
lambda = 50;
robustConfidence=0.9;
% number of observations in each sample
T=100;
% number of portfolios formed by resampling
L=60;
% number of representatives
k=12;
% CVaR scenarios
scenarios=2000;

cvarConfidence=0.95;
% Other cases: [confidence scenarios]
otherConfig=[0.9,10;0.9 1e4;0.99 10;0.99 1e4];
NoCases=size(otherConfig,1);
Casenames=cellstr("CI"+otherConfig(:,1)*100+"Samples"+otherConfig(:,2))';
% TODO
% timesteps, 6 months worth of weeks
periods=26;
%% 3. Construct and rebalance your portfolios
% Here you will estimate your input parameters (exp. returns, cov. matrix, etc) 
% from the Fama-French factor models. You will have to re-estimate your parameters 
% at the start of each rebalance period, and then re-optimize and rebalance your 
% portfolios accordingly.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;

% initialize variables
% portfolio for each method
clear x currentVal NoShares tCost portfValue xT Simulation Buckets ResampleMVOs RobustInfo
shortmethods={'mvo','robust','resampling','diverse','cvar'};
shortmethods=[shortmethods Casenames];
NoMethods=NoMethodsOld+NoCases;
x=cell(1,NoPeriods);
for i=1:NoPeriods
    x{i}=array2table(zeros(n,NoMethods),'RowNames',tickers,'Variablenames',shortmethods);
end
currentVal = zeros(NoPeriods,NoMethods);
NoShares = zeros(n,NoMethods);
tCost = array2table(zeros(NoPeriods-1,NoMethods),'Variablenames',shortmethods);
portfValue = array2table(zeros(size(dates(testStart <= dates),1),NoMethods),'Variablenames',shortmethods);
Simulation = cell(NoCases+1,NoPeriods);
Buckets = cell(1,NoPeriods);
ResampleMVOs=cell(1,NoPeriods);
RobustInfo=cell(1,NoPeriods);
% stores dates range for each period
clear windows SharpesAnte SharpesPost
windows.estimation = cell(1,NoPeriods);
windows.test = cell(1,NoPeriods);
SharpesAnte=zeros(NoMethods,NoPeriods);
SharpesPost=zeros(NoMethods,NoPeriods);
targets=zeros(1,NoPeriods);
Correlations=zeros(1,NoPeriods);

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    currdates = calStart <= dates & dates <= calEnd;
    periodReturns =  returns( currdates, :) ;
    periodFactRet = table2array( factorRet( currdates, :) );
    %     This returns only one row of the current weekly price
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ...
        & dates <= calEnd, :) )';
    % Subset the prices corresponding to the current out-of-sample test
    % period.
    testdates =testStart <= dates & dates <= testEnd;
    periodPrices=adjClose( testdates,:);
    periodPrices = table2array( periodPrices);
    
    windows.estimation{t}=currdates;
    windows.test{t}=testdates;
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        currentVal(1,:) = ones(1,NoMethods)*initialVal;
    else
        currentVal(t,:) = currentPrices'* NoShares;
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % Calculate your initial exp. returns and covariance matrix using the
    % Fama-French factor model. You may find these values as you prefer,
    % the code given below is just an example.
    
    %     Linear Regression
    %     Each column is factors for an asset
    factormodel=FamaFrench(periodReturns,periodFactRet);
    mu=factormodel.mu;
    Q=factormodel.Q;
    %----------------------------------------------------------------------
    
    % Define the target return
    targetRet = mean(mu);
        targets(t)=targetRet;

    x{t}.mvo=MVO(mu, Q, targetRet);
    [x{t}.robust, RobustInfo{i}]=robust(mu, Q, targetRet, periods, lambda,robustConfidence);
    [x{t}.resampling, ResampleInfo]=resampling(mu, Q, targetRet,T,L);
    [x{t}.diverse, BucketInfo]=diverse(mu, Q, targetRet,k);
    [x{t}.cvar , Simulation{1,t}]=cvar(mu, Q, targetRet,scenarios, periodPrices(end,:), periods, cvarConfidence);
    [x{t}.CI90Samples10 , Simulation{2,t}]=...
        cvar(mu, Q, targetRet,otherConfig(1,2), periodPrices(end,:), periods, otherConfig(1,1));
    [x{t}.CI90Samples10000 , Simulation{3,t}]=...
        cvar(mu, Q, targetRet,otherConfig(2,2), periodPrices(end,:), periods, otherConfig(2,1));
    [x{t}.CI99Samples10 , Simulation{4,t}]=...
        cvar(mu, Q, targetRet,otherConfig(3,2), periodPrices(end,:), periods, otherConfig(3,1));
    [x{t}.CI99Samples10000, Simulation{5,t}]=...
        cvar(mu, Q, targetRet,otherConfig(4,2), periodPrices(end,:), periods, otherConfig(4,1));
    
    Buckets{t}=BucketInfo.z;
    Correlations(t)=BucketInfo.totalCorr;
    ResampleMVOs{t}=ResampleInfo.optimals;
    rawx=table2array(x{t})';
    
    %     annualized sharpe ratio
    SharpesAnte(:,t)=rawx*((1+mu).^52-1)./sqrt(diag(rawx*Q*52*rawx'));
    
    % Calculate the optimal number of shares of each stock you should hols
    NoShares = table2array(x{t}).* currentVal(t,:) ./ currentPrices;
    rawportValue=periodPrices * NoShares;
    portfValue(fromDay:toDay,:) = array2table(rawportValue);
    
    realu=(rawportValue(2:end,:)-rawportValue(1:end-1,:))./rawportValue(1:end-1,:);
    realmean=(1+mean(realu)).^52-1;
    realSigma=std(realu)*sqrt(52);
    SharpesPost(:,t)=(realmean./realSigma)';
    
    % Calculate your transaction costs for the current rebalance
    % period. The first period does not have any cost since you are
    % constructing the portfolios for the 1st time.
    
    if t ~= 1
        tCost(t-1,:) = array2table(0.005 * currentPrices'*abs(NoSharesOld - NoShares) );
    end
    
    NoSharesOld = NoShares;
    %------------------------------------------------------------------
    
    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);
end

xT.mvo=cell2mat(cellfun(@(t) t.mvo,x,'UniformOutput', false ));
% [xT.mvolong,xT.mvoshort]=extractLongShort(xT.mvo);
xT.robust=cell2mat(cellfun(@(t) t.robust,x,'UniformOutput', false ));
xT.resampling=cell2mat(cellfun(@(t) t.resampling,x,'UniformOutput', false ));
% [xT.reslong,xT.resshort]=extractLongShort(xT.resampling);
xT.diverse=cell2mat(cellfun(@(t) t.diverse,x,'UniformOutput', false ));
% [xT.divlong,xT.divshort]=extractLongShort(xT.diverse);
xT.cvar=cell2mat(cellfun(@(t) t.cvar,x,'UniformOutput', false ));
xT.additional={cell2mat(cellfun(@(t) t.CI90Samples10,x,'UniformOutput', false )),...
    cell2mat(cellfun(@(t) t.CI90Samples10000,x,'UniformOutput', false )),...
    cell2mat(cellfun(@(t) t.CI99Samples10,x,'UniformOutput', false )),...
    cell2mat(cellfun(@(t) t.CI99Samples10000,x,'UniformOutput', false ))};
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reinitialize testStart
testStart = datetime('2013-01-01');
plotDates = dates(testStart <= dates);

% Reference Portfolio
%   imagine a continousely adjusting market portfolio, always follow market
%   cap weights
portMarket=ones(size(portfValue,1),1)*100;
x_mkt=zeros(n,size(portMarket,1));

currentPrices=table2array( adjClose(plotDates(1) , :) )';
x_mkt(:,1) = mktNoShares .* currentPrices ./ sum( mktNoShares .* currentPrices );
NoShares = x_mkt(:,1) .* portMarket(1) ./ currentPrices;

for thisweek=2:size(plotDates,1)
    %     new prices, portfolio grow
    currentPrices=table2array( adjClose(plotDates(thisweek) , :) )';
    portMarket(thisweek) = currentPrices'* NoShares;
    %     new market weights, rebalance
    x_mkt(:,thisweek) = mktNoShares .* currentPrices ./ sum( mktNoShares .* currentPrices );
    NoShares = x_mkt(:,thisweek) .* portMarket(thisweek) ./ currentPrices;
end
portfValue.market=portMarket;
% Move market to the first one
portfValue=[portfValue(:,end) portfValue(:,1:end-1)];
%% 
% Portfolio Value Plot
%%
names = portfValue.Properties.VariableNames;
rawportfValue=table2array(portfValue);
portValuePlot = figure;
[valuesToPlot,order]=sortSeries(rawportfValue(:,1:NoMethodsOld+1));
line(plotDates, valuesToPlot);
legend(names(order), 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);
drawnow
%% 

AdditionalportValuePlot = figure;
[cvarValues,I]=sortSeries(rawportfValue(:,NoMethodsOld+1:end));
additionalnames=names(NoMethodsOld+1:end);
line(plotDates, cvarValues);
legend(additionalnames(I), 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Value Additional', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);
drawnow
%% 
% Portfolio weights
%%

mktPlot=figure;
plotWeight(x_mkt','Market',tickers);

mvoPlot=figure;
plotWeight(xT.mvo','MVO',tickers);
mvoPlots=figure('units','normalized','outerposition',[0 0 0.5 0.7]);
plotWeightChanges(xT.mvo,n,tickers,'MVO Weights In Each Asset Over Time','Weights')

robustPlot=figure;
plotWeight(xT.robust','Robust',tickers);
robustPlots=figure('units','normalized','outerposition',[0 0 0.5 0.7]);
plotWeightChanges(xT.robust,n,tickers,'Robust Weights In Each Asset Over Time','Weights')

resamplingPlot=figure;
plotWeight(xT.resampling','Resampling',tickers);
resamplingPlots=figure('units','normalized','outerposition',[0 0 0.5 0.7]);
plotWeightChanges(xT.resampling,n,tickers,'Resampling Weights In Each Asset Over Time','Weights')

diversePlot=figure;
plotWeight(xT.diverse','Most-diverse',tickers);
diversePlots=figure('units','normalized','outerposition',[0 0 0.5 0.7]);
plotWeightChanges(xT.diverse,n,tickers,'Diverse Weights In Each Asset Over Time','Weights')

cvarPlot=figure;
plotWeight(xT.cvar','CVaR',tickers);
clear casesFig
casesFig=cell(1,NoCases);
for config=1:NoCases
    casesFig{config}=figure;
    plotWeight(xT.additional{config}',Casenames(config),tickers);
end
drawnow
%% 
% Transaction Cost
%%
transactionCost =tCost;
transactionPlot = figure;
[trans,I] = sortSeries(table2array(transactionCost));
% sort transaction cost by ending value so that legends display in order
line(1:NoPeriods-1, trans)
legend(shortmethods(I), 'Location', 'eastoutside','FontSize',12);
title('Transaction Cost', 'FontSize', 14);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);
drawnow

%% 
% Risk & Return
%%
% 
returns = (rawportfValue(2:end,:)-rawportfValue(1:end-1,:))./rawportfValue(1:end-1,:);
mean_return = (Geomean(1+returns)-1);
stddev= std(returns);
summary=[mean_return;stddev]'*100;
meanstd = array2table(round(summary',3)+"%",...
    'VariableNames',names,'RowNames',{'Mean','Standard Deviation'});
summaryplot=figure;
scatter(summary(:,1),summary(:,2))
text(summary(:,1)*1.01,summary(:,2)*1.01,names)
xlabel('Return %')
ylabel('Standard Deviation %')
drawnow
%% 
% Show how resampling works
%%
resplots=cell(1,NoPeriods);
for t=1:1
    % for t=1:NoPeriods
    resplots{1,t}=figure('units','normalized','outerposition',[0 0 0.5 0.7]);
    subplot(2,1,1)
    plotWeightChanges(ResampleMVOs{t},n,tickers,'Sampling Portfolios For One Period','Weights');
    
    subplot(2,1,2)
    plotWeightChanges(x{t}.resampling,n,tickers,'Result Portfolio','Weights');
end
drawnow

%% 
% Show Representatives for A period
%%
representation=cell(n,NoPeriods);
for t=1:NoPeriods
    z=logical(Buckets{t});
    result=cell(n,1);
    for i=1:n
        result{i}=tickers{z(i,:)};
    end
    representation(:,t)=result;
end
representativeInfo=cell2table(representation,'VariableNames',cellstr("Period"+(1:NoPeriods)),'RowNames',tickers)
writetable(representativeInfo,'representative.csv','WriteRowNames',true);
% get correlation between MVO and other asset
mvocors=corr(returns);
mvocors=mvocors(2,:);
mvocorrPlot=figure;
bar(mvocors(2:6));
xticklabels(names(2:6));
set(gca,'XTickLabelRotation',30);
title("Realized Return Correlation with MVO");
%% 
% Plot Monte Carlo Simulations for Each Period
%%
simulationFigs=cell(1,NoPeriods);
% for t=1:NoPeriods
for t=3:3
    S=Simulation{1,t}.S;
    simulationFigs{t}=figure('units','normalized','outerposition',[0 0 0.8 1]);
    rownum=factor(n);
    rownum=rownum(end);
    for i=1:n
        subplot(rownum,n/rownum,i)
        %         just plot 100 scenarios for performance reason
        line(1:2,reshape(S(:,i,1:100),2,100))
        xticks([1 2])
%         xlabel('Timestep')
        title(tickers(i))
    end
end
drawnow
%% 
% Sharpe Ratio
%%
% Ex Post sharpe ratio comparison across portfolios
ExPostPort=figure;
[valueToPlot,I]=sortSeries(SharpesPost');
plot(valueToPlot);
legend(shortmethods(I),'Location', 'eastoutside')
title('Ex Post Sharpe Ratio')
xlabel('Period')
ylabel('Sharpe')
% Ex Ante
ExAntePort=figure;
[valueToPlot,I]=sortSeries(SharpesAnte');
plot(valueToPlot);
legend(shortmethods(I),'Location', 'eastoutside')
title('Ex Ante Sharpe Ratio')
xlabel('Period')
ylabel('Sharpe')
% Difference
DiffSharpe=figure('units','normalized','outerposition',[0 0 0.7 0.8]);
bar(SharpesAnte-SharpesPost);
xticks(1:NoMethods)
set(gca,'XTickLabelRotation',30);
xlim([0 NoMethods+1])
xticklabels(shortmethods)
title('(Ex Ante - Ex Post) Sharpe')
ylabel('Sharpe')
drawnow
%% 
% CVaR for CVaR portfolio
%%
clear cvmvo cvrobust cvresample cvdiverse cvcvar
cvarPlotsize=[0 0 0.8 0.8];
CVarPlot=figure('units','normalized','outerposition',cvarPlotsize);
[ cvcvar.VaRs,cvcvar.cVaRs ] = getCVaR(Simulation(1,:),xT.cvar, NoPeriods,cvarConfidence );
for t=1:NoPeriods
    assert(abs(cvcvar.VaRs(t)-Simulation{1,t}.VaR)<1e-5,'VaR from optimization should be consistent')
    assert(abs(cvcvar.cVaRs(t)-Simulation{1,t}.cVaR)<1e-5,'CVaR from optimization should be consistent')
end
drawnow
%% 
% CVaR for each porfolio

MVOCvarPlot=figure('units','normalized','outerposition',cvarPlotsize);
[ cvmvo.VaRs,cvmvo.cVaRs ] = getCVaR(Simulation(1,:),xT.mvo, NoPeriods,cvarConfidence);
RobustCvarPlot=figure('units','normalized','outerposition',cvarPlotsize);
[ cvrobust.VaRs,cvrobust.cVaRs ] = getCVaR(Simulation(1,:),xT.robust, NoPeriods,cvarConfidence);
ResampleCvarPlot=figure('units','normalized','outerposition',cvarPlotsize);
[ cvresample.VaRs,cvresample.cVaRs ] = getCVaR(Simulation(1,:),xT.resampling, NoPeriods,cvarConfidence);
DiverseCvarPlot=figure('units','normalized','outerposition',cvarPlotsize);
[ cvdiverse.VaRs,cvdiverse.cVaRs ] = getCVaR(Simulation(1,:),xT.diverse, NoPeriods,cvarConfidence);
drawnow

% CustomizeCvarplots=cell(1,NoCases);
% 
% for i=1:NoCases
%     CustomizeCvarplots{i}=figure('units','normalized','outerposition',[0 0 1 1]);
%     getCVaR(Simulation(1+i,:),xT.additional{i}, NoPeriods,otherConfig(i,1));
% end

%% Compare CVaR for two portfolios
cvComparePlot=figure('units','normalized','outerposition',cvarPlotsize);
plotPairCVaRs({Simulation(1,:),Simulation(1,:)},{xT.cvar,xT.mvo},...
    shortmethods([5 1]),6,[cvarConfidence cvarConfidence]);
cvCIComparePlot=figure('units','normalized','outerposition',cvarPlotsize);
plotPairCVaRs({Simulation(3,:),Simulation(5,:)},{xT.additional{2},xT.additional{4}},...
    Casenames([2 4]),6,[otherConfig(2,1) otherConfig(4,1)]);
%% 
% In some cases CVaR optimization may not produce the least CVaR, this is 
% because other portfolios allow short-sell. In general, CVaR optimization obtains 
% lowest CVaR.

% Compare across portfolios
[CVaRs,I]=sortSeries([cvmvo.cVaRs' cvrobust.cVaRs' cvresample.cVaRs' cvdiverse.cVaRs' cvcvar.cVaRs']);
% VaRs=[cvmvo.VaRs' cvrobust.VaRs' cvresample.VaRs' cvdiverse.VaRs' cvcvar.VaRs'];
CVaRAllPlots=figure('units','normalized','outerposition',[0 0 0.6 0.7]);

plot(CVaRs);
xlabel('Periods')
ylabel('CVaR')
legend(shortmethods(I),'Location','eastoutside')
drawnow
%%
% additional report assisting plots
pricesPlot=plotTimeTable(adjClose);
drawnow
%% 
% Store All Plots
%%

figurestoprint={pricesPlot,portValuePlot,AdditionalportValuePlot,...
    mktPlot,mvoPlot,mvoPlots,robustPlot,robustPlots,resamplingPlot,resamplingPlots,...
    diversePlot,diversePlots,cvarPlot,transactionPlot,summaryplot,resplots{1},...
    simulationFigs{3},ExPostPort,ExAntePort,DiffSharpe,...
    CVarPlot,cvComparePlot,CVaRAllPlots,cvCIComparePlot,mvocorrPlot };
filenames={'AllStockPrice','PortfolioVal','additionalPortfolioVal','market','MVO','MVObar'...
'robust','robustbar','resample','resamplebar','diverse','diversebar','cvar','transaction',...
'overview','howResampling','howMonteCarlo','ExPost','ExAnte',...
'sharpeDiff','cvCVaR','cvCompare','CVaRAll','CV9099Compare','MVOCorrelation'};
fprintf("Printing Graphs Please wait...\n")
for i=1:size(figurestoprint,2)
    print(figurestoprint{i},filenames{i},'-dpng','-r250');
end
for c=1:NoCases
    print(casesFig{c},Casenames{c},'-dpng','-r250');
end
fprintf("Done Printing Graphs\n")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End