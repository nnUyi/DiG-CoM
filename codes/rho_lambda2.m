function [Rho_and_lambda2] = rho_lambda2()
%%       Get Rho and lambda2 for group images
rain_dir = 'E:\datasets\rain_data_train_Light\rain\';
title = 'Rain100L-train';
fm = 'norain-*x2.png'; % format
files = dir([rain_dir '\' fm]);
img_num = length(files);

Rho_and_lambda2 = zeros(900, 2); %% rho lambda2
for i=1:900
    img_name=[rain_dir sprintf('norain-%dx2.png', i)]
    % img_name = [rain_dir files(i).name]
    rainy_img = double(imread(img_name))./255;
    [Y, ~, ~] = rgb2yuv(rainy_img);
    coarse_mask = RainDetectConv(rainy_img, 9, 7, 0.01, 0.12);
    [~, ~, ~, ~, propotion, ~, ~] = MaskModifyImpro(coarse_mask, Y, 5, 8, 3, 0, 0);
    lambda2 = min(0.01, 0.0045*(3.65)/(propotion*100)^2);
    Rho_and_lambda2(i, 1) = propotion;
    Rho_and_lambda2(i, 2) = lambda2;
end
save([rain_dir 'Rho_lambda2_' title '.mat'], 'Rho_and_lambda2');
end