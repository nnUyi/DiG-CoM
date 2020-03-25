function [Mask, theta, Theta_cluster, amount, propotion, M1, M2] = MaskModifyImpro(Mark, Y,  T1, T2, T3, averTheta, flag)
%  Use the properties of rain to modify mark got by first step.
%  Line like shape and angle
%  amount represents the amount of rain
%  averTheta: average angle of whold image
%  flag: whether to use averTheta
%  IF YOU WANT TO CHANGE THE ANGKE MANUALLY, CHANGE THE MOST_FREQUENT_ANGLE
%  

%%                Parameters
  % Mark: input binary image, 0 for rain pixels and 1 for non-rain pixels
  % T1: threshold parameter one, for line like shape detect
  % T2: threshold parameter two, for angle detect
  % T3: threshold parameter three, for length of rain streak

tic;
Mask = 1 - Mark;  % 1 for foreground pixels as rain pixels
amount = 0.0;
% imshow(Mask);
%title('Initial Mark');
% pause;
[H, W] = size(Mask);
M1 = zeros(H, W);
M2 = zeros(H, W);
%%                1. Get connected components from Mark
[L, num] = bwlabel(Mask, 8);  % use 8-adjacency
fprintf('There are %d connected components to handle in total.\n', num);
Theta_cluster = 180.*ones(num, 1);
cnt = 0;
Lambda = zeros(num, 1);
%%               2. Testing all connected components
EulerTopo = regionprops(L, 'EulerNumber');   %% No holes in rain streaks
PixelIdxLists = regionprops(L, 'PixelIdxList');
PixelLists = regionprops(L, 'PixelList');
for iter = 1:num
    Indices = PixelIdxLists(iter).PixelIdxList;
    if size(Indices) == 1
        continue
    end
    if (EulerTopo(iter).EulerNumber) ~= 1 %% hole in rain streaks
        continue
    end
    location = PixelLists(iter).PixelList;
    location = location(:, [2 1]); %????
    location = location - mean(location, 1);
    [V, D] = eig(location'*location);  % D is diagonal matrix, V is the eigen matrix
    lambda1 = abs(D(1, 1));
    lambda2 = abs(D(2, 2));
    angle = atand(V(2, 2)/V(1, 2));
    if (lambda2/(1e-5+lambda1) > T1) && (lambda1 > T3) && (abs(angle)<=45 && lambda2<1e6 && lambda1 < 1e3)
        Theta_cluster(iter, 1) = angle;
        Lambda(iter, 1) = lambda1;
        % Mask(Indices) = 1;
    else
        Mark(Indices) = 0;
    end
    if lambda2 > 1e6
        M1(Indices) = 1;
    end
end
if flag == 1
    most_frequent_angle = averTheta;
else
    [Counts, Value] = hist(Theta_cluster(Theta_cluster~=180), ceil(num/8));  % every angle 8 rain streaks in average
    Counts = conv1d(Counts, ones(1, 11)/11); % use 1d convolution to get the precise angle, cause rain streaks take the same angle, thus the neighbour angle value should be close
    % 卷积提取角度,雨线角度接近,因此卷积后,会把角度统计值放大,避免选错角度
    [~, idx] = max(Counts);
    most_frequent_angle = min(180, Value(idx)); 
end
%%                     3. Mark modify
Theta = zeros(num, 1);
for iter = 1:num
    Indices = PixelIdxLists(iter).PixelIdxList;
    angle = Theta_cluster(iter);
    if (sum(sum(Indices)) < 3)
        Mask(Indices) = 0;
        % M2(Indices) = 1;
        continue;
    end
    if (abs(angle-most_frequent_angle) > T2)
        Mask(Indices) = 0; % non_rain pixels
    
    else
        amount = amount + Lambda(iter, 1);
        cnt = cnt + 1;
        Theta(cnt) = angle;
    end
    if angle ~= 180 && abs(angle-most_frequent_angle) > T2
        M2(Indices) = 1;
    end
end
        
Theta_cluster = Theta(1:cnt, 1);
[Counts, Value] = hist(Theta_cluster, ceil(cnt/8));  % every angle 8 rain streaks in average 
% plot(Value, Counts, 'r');
[~, idx] = max(Counts);
% pause;
theta = Value(max(1, idx-1):min(idx+1, size(Value, 2)));
weight = Counts(max(1, idx-1):min(idx+1, size(Counts, 2)));
if sum(weight) == 0  % no rain streaks, set theta to 180 degree
    theta = 180;
    amount = 0.0;
    fprintf('No rain streaks detected, set theta to 180 degree.\n');
else
    theta = sum(weight.*theta)/sum(weight);
    % theta = sum(Counts.*Value)/sum(Counts);
end
fprintf('The average angle is %.4f.\n', theta);
% imshow(Mask);
%%propotion = sum(sum(Mask.*Y))/(size(Mask, 1)*size(Mask, 2));
propotion = sum(sum(Mask.*Y))/(size(Mask, 1)*size(Mask, 2));
% title('Final Mark');
% pause;
Mask = 1 - Mask;
t2 = toc;
fprintf('markmodify using regionprops time consuming: %.6f\n', t2);
end