function [ x_optimal, otherInfo ] = resampling( mu, Q, targetRet, T, L )
% T: sample sizes
% L: number of generated portfolios
% otherInfo.u constains L elements, each is the u of a sample
% otherInfo.Q constains L elements, each is the Q of a sample
% otherInfo.x constains n x L elements, each column is a portfolio generated
% from the particular sample population
n = size(Q,2);

Qsample = cell(L, 1);
meanM = zeros(n, L);
optimal = zeros(n, L);

for i = 1:L
    RandRet = mvnrnd(mu, Q, T);
    meanM(:,i) = (Geomean(1+RandRet)-1)';
    Qsample{i} = cov(RandRet);
    optimal(:,i) = MVO(meanM(:,i),Qsample{i},targetRet); 
end
x_optimal = mean(optimal, 2);
otherInfo.u=meanM;
otherInfo.Q=Qsample;
otherInfo.optimals=optimal;
end