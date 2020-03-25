function [cols, Mark, thetas, propotions] = preProcess(img)
%%    decide how many rows and cols to divide
%%    initial is 2 by 2 pieces
%%    return number pieces and Mask

thetas = zeros(1, 2);
propotions = zeros(1, 2);
W = size(img, 2);
halfW = floor(W/2);
cols = 1;

SE = strel('arbitrary', ones(3, 3));
[Y, ~, ~] = rgb2yuv(img);
  % YUV = rgb2ycbcr(I);
  % Y = YUV(:, :, 1);
  Mark = RainDetectConv(img, 9, 7, 0.01, 0.12);
  T2 = 25;
  MarkL = Mark(:, 1:halfW+1);
  MarkR = Mark(:, halfW:end);  %% ÓÐ½»²æ
  % [MarkL, ThetaL, ~, ~, propotionL] = MarkModify(MarkL, Y(:, 1:halfW+1), 5, T2, 3, 0, 0);
  [MarkL, ThetaL, ~, ~, propotionL] = MaskModifyImpro(MarkL, Y(:, 1:halfW+1), 5, T2, 3, 0, 0);
  % [MarkR, ThetaR, ~, ~, propotionR] = MarkModify(MarkR, Y(:, halfW:end), 5, T2, 3, 0, 0);
  [MarkR, ThetaR, ~, ~, propotionR] = MaskModifyImpro(MarkR, Y(:, halfW:end), 5, T2, 3, 0, 0);
  
  if abs(ThetaL - ThetaR) > 8 || ThetaL*ThetaR < 0
      cols = 2;
      thetas(1) = ThetaL;
      thetas(2) = ThetaR;
      propotions(1) = propotionL;
      propotions(2) = propotionR;
  end
  if propotionL/propotionR < 0.3 || propotionL/propotionR > 3
      cols = 2;
      propotions(1) = propotionL;
      propotions(2) = propotionR;
      thetas(1) = ThetaL;
      thetas(2) = ThetaR;
  end
  thetas(1) = (ThetaL+ThetaR)/2;
  propotions(1) = (propotionL+propotionR)/2;
  Mark(:, 1:halfW) = MarkL(:, 1:halfW);
  Mark(:, halfW+1:end) = MarkR(:, 2:end);
  Mark = 1-imclose(1-Mark, SE);
  figure
  imshow(Mark);
  title('Mark');
end 