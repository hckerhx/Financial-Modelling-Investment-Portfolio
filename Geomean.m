function mean = Geomean(x)
if ispc
    %     on windows machines this is geo_mean
    mean=geo_mean(x);
else
    %     on mac its geomean
    mean=geomean(x);
end
end