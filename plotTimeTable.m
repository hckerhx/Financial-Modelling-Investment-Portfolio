function [fig] = plotTimeTable(data)
% plot array, each column is a series
fig=figure;
x = data.Properties.RowTimes;
legends = data.Properties.VariableNames;
[Values,I]=sortSeries(table2array(data));
line(x,Values);
legend(legends(I),'Location', 'eastoutside')
end