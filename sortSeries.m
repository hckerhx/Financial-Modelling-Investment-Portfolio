function [ sortedseries,I ] = sortSeries( series )
% sort serieses by the end value
[~, I]=sort(-series(end,:));
sortedseries=series(:,I);
end

