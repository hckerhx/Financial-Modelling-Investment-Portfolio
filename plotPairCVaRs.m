function [] = plotPairCVaRs(Simulations,xs, names, NoPeriods,confidences)
% plot pair CVars on one graph
Simulation1=Simulations{1};
Simulation2=Simulations{2};
[scenarios1,n]=size(Simulation1{1}.r);
[scenarios2,~]=size(Simulation2{1}.r);

% f(x,ys)
fxys1=zeros(scenarios1,NoPeriods);
fxys2=zeros(scenarios2,NoPeriods);

threshold1=uint16(scenarios1*confidences(1));
threshold2=uint16(scenarios2*confidences(2));

numrows = factor(NoPeriods);
numrows=numrows(end);
numcols=NoPeriods/numrows;

% bin width
BW=0.8/50;
for t=1:NoPeriods
    r1=Simulation1{t}.r;
    r2=Simulation2{t}.r;

    fxys1(:,t)=-r1*xs{1}(:,t);
    fxys2(:,t)=-r2*xs{2}(:,t);

    subplot(numrows,numcols,t)
    
    histogram(fxys1(:,t),'BinWidth',BW)
    hold on
    histogram(fxys2(:,t),'BinWidth',BW)

    xlabel('Loss')
    xlim([-0.4 0.4])
    title("Profit and loss Period "+t)
    sortedScenarios1=sort(fxys1(:,t));
    sortedScenarios2=sort(fxys2(:,t));

    VaRs1=sortedScenarios1(threshold1);
    cVaRs1=mean(sortedScenarios1(threshold1+1:end));
    
    VaRs2=sortedScenarios2(threshold2);
    cVaRs2=mean(sortedScenarios2(threshold2+1:end));
    linescalor=0.7/sqrt(n);
    line([VaRs1, VaRs1], [0, scenarios1*linescalor],'Color','r')
    line([VaRs2, VaRs2], [0, scenarios2*linescalor],'Color','g')
    line([cVaRs1, cVaRs1], [0, scenarios1*linescalor],'Color','b')
    line([cVaRs2, cVaRs2], [0, scenarios2*linescalor],'Color','black')
    legend([names,"VaR "+names,"CVaR "+names],'FontSize',12)
end    
    
end

