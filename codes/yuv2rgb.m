function [res] = yuv2rgb(Input)
[c] = size(Input, 3);
res = zeros(size(Input));
if c == 1
    Input = repmat(Input, [1 1 3]);
end
Y = Input(:, :, 1);
U = Input(:, :, 2);
V = Input(:, :, 3);

R = Y+1.403.*V;
G = Y-0.344.*U-0.714.*V;
B = Y+1.773.*U;
res(:, :, 1)=R;
res(:, :, 2)=G;
res(:, :, 3)=B;
end