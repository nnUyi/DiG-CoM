% This function is used to get the average theta of the whole image, for
% the sake of the initial theta

function [averTheta] = GetTheta(I, T1, T2)
Mark = RainDetectConv(I, 7, 7, 0.01, 0.12);
Mark = 1 - Mark;
[L, num] = bwlabel(Mark, 8);
Theta_cluster = 180.*ones(num, 1);
for iter = 1:num
    Indices = find(L == iter);
    location = zeros(size(Indices, 2));  % location of all pixels in the same connected components
    for i = 1:size(Indices)
        [r, c] = ind2sub(size(Mark), Indices(i));
        location(i, 1) = r;
        location(i, 2) = c;
    end
    location = location - mean(location);
    [V, D] = eig(location'*location);  % D is diagonal matrix, V is the eigen matrix
    lambda1 = abs(D(1, 1));
    lambda2 = abs(D(2, 2));
    angle = atand(V(2, 2)/V(1, 2));
    if (lambda2/(1e-5+lambda1) > T1) && (lambda1 > T2) && (abs(angle)<=45)
        Theta_cluster(iter, 1) = angle;
    else
        Mark(Indices) = 0;
    end
end
[Counts, Value] = hist(Theta_cluster(Theta_cluster~=180), ceil(num/8));
averTheta = sum(Counts.*Value)./sum(Counts);

end