function [] = plotWeightChanges(x,n,xlabels,ctitle,cylabel )
bar(x)
xticks(1:n)
xlim([1-1 n+1])
xticklabels(xlabels)
title(ctitle)
ylabel(cylabel)
end

