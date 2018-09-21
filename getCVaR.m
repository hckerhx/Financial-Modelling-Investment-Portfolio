function [ VaRs,cVaRs ] = getCVaR(Simulation,x, NoPeriods,confidence )
% Simulation: cell with each element being one simulation struct,
% r:returns, S: prices
% x: n x NoPeriods, portfolio
%
[scenarios,n]=size(Simulation{1}.r);
% f(x,ys)
fxys=zeros(scenarios,NoPeriods);
VaRs=zeros(1,NoPeriods);
cVaRs=zeros(1,NoPeriods);
threshold=uint16(scenarios*confidence);
% for subplotting, defind rows and cols
numrows = factor(NoPeriods);
numrows=numrows(end);
numcols=NoPeriods/numrows;
for t=1:NoPeriods
    r=Simulation{t}.r;
    fxys(:,t)=-r*x(:,t);
    subplot(numrows,numcols,t)
    histogram(fxys(:,t))
    xlabel('Loss')
    xlim([-0.4 0.4])
    title("Profit and loss Period "+t)
    sortedScenarios=sort(fxys(:,t));
    %     VaR is simply the 95th largest loss
    VaRs(t)=sortedScenarios(threshold);
    %     cVaR is the expectation (mean) of all losses above the threshold
    cVaRs(t)=mean(sortedScenarios(threshold+1:end));
    
    line([VaRs(t), VaRs(t)], [0, scenarios*0.7/sqrt(n)],'Color','r')
    line([cVaRs(t), cVaRs(t)], [0, scenarios*0.7/sqrt(n)],'Color','b')
    legend({'Loss','VaR','CVaR'},'FontSize',12)
end

end

