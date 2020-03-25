clc;
clear all;
close all;
save_path = 'F:\Education\temporal_store';
mode='synthetic';
% quality_res_ours = zeros(100, 3);
% addpath('F:\Education\deeplearningPractice\Derain\Real_Rain\all_rain');
%addpath('F:\Education\deeplearningPractice\Derain\Rain100H');
addpath('E:\datasets\rain12');
% addpath('E:\datasets\Rain100L');
% addpath('F:\Education\deeplearningPractice\Derain\real_test_1000\rain');
% addpath('F:\Education\deeplearningPractice\Derain\real_test_1000\gt');
% addpath('F:\Education\deeplearningPractice\Derain\real_world_rainy_dataset\real_world_rainy_dataset\rain'); % natural rain
% addpath('E:\datasets\Real_Internet\Real_Internet');
addpath 'F:\Education\temporal_store';
Index_Set = IndexExtract('E:\datasets\rain12\');
fprintf('There are %d images to process.\n', Index_Set);
img_num = Index_Set;
addpath(save_path);
if ~exist('quality_result_ours.mat', 'file')
    quality_res_ours = zeros(img_num, 3);
else
    load('quality_result_ours.mat');
end
for index = 1:img_num
    close all;
    if exist(['00' num2str(index) '_Background_ours.png'], 'file')
        continue;
    else
        if index < 10
            img_name = ['00' num2str(index) '_in.png'];  %rain12
            gt_name = ['00' num2str(index) '_GT.png'];   %rain12
            % img_name = ['rain-00' num2str(index) '.png'];  %rain100L
            % gt_name = ['norain-00' num2str(index) '.png'];
            % img_name = ['00' num2str(index) '.jpg'];  %natural rain
            %img_name = ['00' num2str(index-1) '.png']; % real test 1000
           % gt_name = ['00' num2str(index-1) 'gt.png']; % real test 1000
           % img_name = ['rain-00' num2str(index) '.png']; % real internet
        elseif index < 100
            img_name = ['0' num2str(index) '_in.png'];    %rain12
            gt_name = ['0' num2str(index) '_GT.png'];     %rain12
            % img_name = ['rain-0' num2str(index) '.png'];    %rain100L
            % gt_name = ['norain-0' num2str(index) '.png'];   %rain100L
            % if index==14
            %     img_name = ['0' num2str(index) '.png'];    %natural rain
            % else
            %      img_name = ['0' num2str(index) '.jpg'];    %natural rain
            % end
            %img_name = ['0' num2str(index-1) '.png']; % real test 1000
           % gt_name = ['0' num2str(index-1) 'gt.png']; % real test 1000
            % img_name = ['0' num2str(index) '.jpg']; % natural rain images
            %img_name = ['rain-0' num2str(index) '.png']; % real internet
        else
            % img_name = ['rain-' num2str(index) '.png'];    %rain100L
            % gt_name = ['norain-' num2str(index) '.png'];   %rain100L
           % img_name = [num2str(index-1) '.png']; % real test 1000
           % gt_name = [num2str(index-1) 'gt.png']; % real test 1000
            % img_name = [num2str(index) '.jpg']; % natural rain images
            % img_name = ['rain-' num2str(index) '.png']; % real internet
        end
        I = double(imread(img_name))./255;
        BT = double(imread(gt_name))./255;
        figure;
        imshow(I);
        size_of_blocks = [1 1];
        % size_of_blocks = input('Please input the size of blocks according to the raw image, examples like: [1 1] or [1 2]: ');
        if isempty(size_of_blocks)
            size_of_blocks = [1 1];
        end
        fprintf('**************** Program begin ***************\n');
        tic;
        [YB, YR, B, R, ~, ~, ~] = DerainByWindow(I, size_of_blocks(1), size_of_blocks(2), mode);
        t2 = toc;
        quality_res_ours(index, 1) = ssim(B, BT);
        quality_res_ours(index, 2) = psnr(B, BT);
        quality_res_ours(index, 3) = t2;
        %quality_res_ours(index, 1) = niqe(B);
        %quality_res_ours(index, 2) = t2;
        imwrite(YB./255, [save_path '\00' num2str(index) '_YB.png']);
        imwrite(B, [save_path '\00' num2str(index) '_Background_ours.png']);
        imwrite(YR+0.5, [save_path '\00' num2str(index) '_Rain_ours.png']);
        fprintf('---------------- The %d image has finished in (%.3fs) -------------------\n', index, t2);
        fprintf('\n\n');
        save([save_path '\quality_result_ours.mat'], 'quality_res_ours');
    end
end
save([save_path '\quality_result_ours_1e-4_1e-4.mat'], 'quality_res_ours');