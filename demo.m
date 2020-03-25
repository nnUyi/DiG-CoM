img_format = 'rain-*.png';
img_path = './imgs/';
mode = 'real'; % synthetic
addpath './codes';
dirOut = dir(fullfile(img_path, img_format));  
Names = {dirOut.name};
len = size(Names, 2);
for i=1:len
    img_name = [img_path Names{i}];
    fprintf('image name: %s\n', Names{i});
    img = double(imread(img_name))./255;
    [~, ~, B, ~, ~, ~] = DerainByWindow(img, 1, 1, mode);
    figure;
    imshow(img);
    title('rainy image');
    figure;
    imshow(B);
    title('clean background');
    fprintf('press Enter to continue....\n');
    imwrite(B, ['./results/' 'clean-no' Names{i}]);
    pause;
end
