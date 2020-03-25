function [res] = conv1d(a, ker)
%% One D convolution operation
if numel(a) < numel(ker)
    res = a;
else
    m = numel(a);
    f = numel(ker);  
    p = (f-1)/2;
    tmpa = zeros(1, m+f-1);
    tmpa(1, 1:p) = a(1, 1);
    tmpa(1, p+1: m+p) = a;
    tmpa(m+p+1: m+f-1) = a(1, end);
    res = zeros(1, m);
    for i=1:f
        res = res + ker(1, i).*tmpa(1, i: m-1+i);
    end
end
end