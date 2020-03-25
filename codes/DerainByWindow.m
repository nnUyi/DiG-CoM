function [YB, YR, B, R, E, B_initial, R_initial, Mark] = DerainByWindow(Input, row_num, col_num, mode)
%%                 Parameters
  % derain by dividing picture to row_num*col_num blocks(patches)
  % row_num: row blocks number
  % col_num: col blocks number
  % IF YOU WANT TO CHANGE LAMBDA2, JUST CHANGE IT MANUALLY.

if size(Input, 3) == 1
    Input = repmat(Input, [1 1 3]);
end
[col_num, ~, ~, ~] = preProcess(Input);
fprintf('columns is : %d\n', col_num);
if (row_num==1) && (col_num==1)
    flag = 0;
    averTheta = 0;
else
    flag = 1;
    averTheta = GetTheta(Input, 3, 3);
end

bias_lambda2 = 0;

if strcmp(mode, 'real')
    bias_lambda2 = 2;
end

[Y, ~, ~] = rgb2yuv(Input);
YB = zeros(size(Y));
YR = zeros(size(Y));
[m, n, ~] = size(Input);
if m > 1000
    row_num = 2;
end
if n > 1000
    col_num = 2;
end
fprintf('cut image into pieces of shape (%d, %d)', row_num, col_num);
B = zeros(size(Input));
E = zeros(size(Input, 1), size(Input, 2));
B_initial = zeros(size(Input));
R = zeros(size(Input));
R_initial = zeros(size(Input));
p_row = floor(m/row_num);  % row length of each block
p_col = floor(n/col_num);  % column length of each block
Mark = zeros(size(Y));
for i=1: row_num-1
    for j=1:col_num-1
        I = Input(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col, :);  % (i, j) block
        [tmpYB, tmpYR, tmpB, tmpR, tmpE, Bappro, Rappro, tmpMark] = DerainBlock(I, averTheta, flag, bias_lambda2); 
        YB(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col) = tmpYB;
        YR(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col) = tmpYR;
        B(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col, :) = tmpB;
        R(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col, :) = tmpR;
        B_initial(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col, :) = Bappro;
        R_initial(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col, :) = Rappro;
        E(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col) = tmpE;
        Mark(1+(i-1)*p_row:i*p_row, 1+(j-1)*p_col:j*p_col) = tmpMark;
    end
end

for j=1:col_num-1
    I = Input(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col, :);  % (row_num, j) block
    [tmpYB, tmpYR, tmpB, tmpR, tmpE, Bappro, Rappro, tmpMark] = DerainBlock(I, averTheta, flag, bias_lambda2); 
    YB(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col) = tmpYB;
    YR(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col) = tmpYR;
    B(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col, :) = tmpB;
    R(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col, :) = tmpR;
    B_initial(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col, :) = Bappro;
    R_initial(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col, :) = Rappro;
    E(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col) = tmpE;
    Mark(1+(row_num-1)*p_row:m, 1+(j-1)*p_col:j*p_col) = tmpMark;
end

for i=1:row_num-1                                
    I = Input(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n, :);  % (i, col_num) block
    [tmpYB, tmpYR, tmpB, tmpR, tmpE, Bappro, Rappro, tmpMark] = DerainBlock(I, averTheta, flag, bias_lambda2); 
    YB(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n) = tmpYB;
    YR(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n) = tmpYR;
    B(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n, :) = tmpB;
    R(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n, :) = tmpR;
    B_initial(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n, :) = Bappro;
    R_initial(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n, :) = Rappro;
    E(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n) = tmpE;
    Mark(1+(i-1)*p_row:i*p_row, 1+(col_num-1)*p_col:n) = tmpMark;
end

I = Input(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n, :);
[tmpYB, tmpYR, tmpB, tmpR, tmpE, Bappro, Rappro, tmpMark] = DerainBlock(I, averTheta, flag, bias_lambda2); 
YB(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n) = tmpYB;
YR(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n) = tmpYR;
B(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n, :) = tmpB;
R(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n, :) = tmpR;
B_initial(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n, :) = Bappro;
R_initial(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n, :) = Rappro;
E(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n) = tmpE;
Mark(1+(row_num-1)*p_row:m, 1+(col_num-1)*p_col:n) = tmpMark;
end


function [YB, YR, B, R, E, Bappro, Rappro, Mark, M1, M2] = DerainBlock(I, averTheta, flag, bias_lambda2)
%%             Paramters
  % derain a block with input I
  % I is the input block image
  % YB is the Y value of final background
  % B is final background
  % R is the final rain layer
  % Bappro is the initial background got by bilateral filter
  % Rappro is the initial rain layer got by bilateral filter
  SE = strel('arbitrary', ones(3, 3));
  [Y, U, V] = rgb2yuv(I);
  %im = rgb2ycbcr(I);
  %Y = im(:, :, 1);
  %%U = im(:, :, 2);
  %V = im(:, :, 3);
  % YUV = rgb2ycbcr(I);
  % Y = YUV(:, :, 1);
  fliplr_flag = false;
  Mark = RainDetectConv(I, 9, 7, 0.01, 0.5);
  if flag == 1
      T2 = 25;   % not noly one block, range from -25 to 25
  else
      T2 = 8;
  end
  % [Mark, Theta, ~, amount, propotion, M1, M2] = MarkModify(Mark, Y, 5, T2, 3, averTheta, flag);
  [Mark, Theta, ~, amount, propotion, M1, M2] = MaskModifyImpro(Mark, Y, 5, T2, 3, averTheta, flag);
  Mark = 1-imclose(1-Mark, SE);
  imshow(Mark)
  title(Mark);
  if Theta == 180   % no rain streaks situation
      YB = Y;
      YR = zeros(size(Y));
      B = I;
      R = zeros(size(I));
      Bappro = B;
      Rappro = R;
      E = zeros(size(Y));
  else
      if Theta > 0  % theta bigger than 0, fliplr the image; theta bigger than 0, means north-west to south-east
          I = fliplr(I);
          Mark = fliplr(Mark);
          Theta = -1*Theta;
          Y = fliplr(Y);
          fliplr_flag = true;
      end
      fprintf('The amount of rain in this block is: %.5f.\n', amount);
      fprintf('The propotion of rain in this block is: %.4f%%.\n', propotion*100);
      lambda2 = min(0.01, 0.0045*(3.65)/(propotion*100+bias_lambda2)^2);
      % lambda2 = min(0.035, 0.015*4.4/(propotion*100)^2);
      fprintf('the parameter lambda2 is %.4f.\n', lambda2);
      % F = BilateralMark(I, 15, Mark, 500, 8.7, 1);
      F = MeanMask(I, 21, 21, Mark);
      Rappro = (I-F).*(1-Mark);
      Bappro = I - Rappro;
      % Bappro = I.*Mark + F.*(1-Mark);
      [YR, ~, ~] = rgb2yuv(Rappro);
      % K = 1/9.*[1 1 1; 1 1 1; 1 1 1];
      % YB = imfilter(Y-YR, K, 'symmetric', 'conv');
      % YE = Y-YR-YB;
      [~, w1, w2] = RotateFilter(-Theta*pi/180);
      K = [0 0 0; 0 1 0; -w2 -w1 0]*0.5;
      % YE = imfilter(Y, K, 'symmetric', 'conv')/2;
      % YB = imfilter(Y-YR, 1/9*ones(3, 3), 'symmetric', 'conv');
      YE = imfilter(Y, K, 'symmetric', 'conv');
      K2 = [0 0 0; 0 1 0; w2 w1 0]*0.5;
      YB = imfilter(Y-YR, K2, 'symmetric', 'conv');
      YR = imfilter(YR, K2, 'symmetric', 'conv');
      % YUVR = rgb2ycbcr(Rappro);
      % YR = YUVR(:, :, 1);
      % Weights = (1-M2) + M2;
      % [B, R, E] = DerainSPBregmanEdgeL1(Y.*255, (Y-YR).*255, YR.*255, 3, lambda2,10, 0.5, 1.0, w1, w2, 1e-4);
      [B, R, E] = derainWithEdgeTrue(Y.*255, YB.*255, YR.*255,1.45, lambda2, 4.1, 0.005, 0.8, w1, w2, 3e-3, YE.*255);
      % [B, R, E] = derainWithEdge(Y.*255, (Y-YR).*255, YR.*255, 0.72, 1e-4, 1.9, 0.5, 0.5, w1, w2, 1e-4);
      % 0.72, lambda2, 1.9, 0.5, 0.5
      % 0.6 0.001 1.5 0.3 0.2
      % 0.8       1.75 0.5 0.6
      % 0.65      1.75 0.3 0.4
      % [B, R] = DerainSPBregman(Y.*255, (Y-YR).*255, YR.*255, 3, lambda2, 9.1, w1, w2, 1e-3);
      % 3 0.015 10.1
      % E = zeros(size(Y));
      if fliplr_flag 
          I = fliplr(I);
          Bappro = fliplr(Bappro);
          Rappro = fliplr(Rappro);
          B = fliplr(B);
          R = fliplr(R);
          E = fliplr(E);
          Y = fliplr(Y);
          Mark = fliplr(Mark);
      end
      R = R./255;
      R = max(0, R);
      R = min(Y, R);
      YB = B;
      YR = R;
     % YR = (1-M1).*YR;
      P = zeros(size(I));
      P(:, :, 2) = U;
      P(:, :, 3) = V;
      % P(:, :, 1) = (YB+E)./255;
      P(:, :, 1) = (Y-YR);
      % YUV(:, :, 1) = (YB+E)./255;
      B = yuv2rgb(P);
      % B = ycbcr2rgb(P);
      % B = ycbcr2rgb(YUV);
      P = zeros(size(I)).*1e-10;
      P(:, :, 1) = R;
      % YUV(:, :, 1) = R./255;
      R = yuv2rgb(P);
      % R = ycbcr2rgb(P);
      % R = ycbcr2rgb(YUV);
      close all;
      %figure;
      %imshow(Y);
      %title('Y channel of source rain image');
      %figure;
      %imshow(YB./255);
      %title('smooth part of background');
      %figure;
      %imshow(Y-YR);
      %title('final result');
  end
end
