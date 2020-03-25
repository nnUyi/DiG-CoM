function [Y, U, V] = rgb2yuv(Input)
[c] = size(Input, 3);
if c == 1
    Input = repmat(Input, [1 1 3]);
end
R = Input(:, :, 1);
G = Input(:, :, 2);
B = Input(:, :, 3);

Y = 0.299.*R + 0.587.*G + 0.114.*B;
U = - 0.169.*R - 0.331.*G + 0.5.*B;
V = 0.5.*R - 0.419.*G - 0.081.*B;
end