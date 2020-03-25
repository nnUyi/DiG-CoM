function [len] = IndexExtract(img_path)

dirOut = dir(fullfile(img_path, '*.png'));  % 获取文件夹下文件名
Names = {dirOut.name};  % 获取文件名的cell
len = size(Names, 2);
end