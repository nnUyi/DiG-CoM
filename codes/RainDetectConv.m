function [Mark] = RainDetectConv(I, w1, w2, mu, epsilon)
%%                   Parameters
  % I: input image, contains 3 channels: R, G, B
  % w1: size of window, length, window change to rectangle
  % w2: size of window, height, window change to rectangle
  % mu: empirical value
  % epsilon: u-v space distance criterion

%%                   Derive Kernel 
  %%                 Initialize Kernel
  K_Left_Up = zeros(2*w1-1, 2*w2-1);
  K_Right_Up = zeros(2*w1-1, 2*w2-1);
  K_Center = 1/(w1*w2).*ones(w1, w2);
  K_Left_Down = zeros(2*w1-1, 2*w2-1);
  K_Right_Down = zeros(2*w1-1, 2*w2-1);
  %%                 Construct Kernel for 5 windows
  K_Left_Up(1:w1, 1:w2) = 1/(w1*w2).*ones(w1, w2);
  K_Right_Up(1:w1, w2:2*w2-1) = 1/(w1*w2).*ones(w1, w2);
  K_Left_Down(w1:2*w1-1, 1:w2) = 1/(w1*w2).*ones(w1, w2);
  K_Right_Down(w1:2*w1-1, w2:2*w2-1) = 1/(w1*w2).*ones(w1, w2);
  
%%                   Give Mark of rain pixels
  % zero for rain pixels
  [m, n, c] = size(I);
  if c==1
      I = repmat(I, [1 1 3]);
  end
  Mark = zeros(m, n);
  I_LU = imfilter(I, K_Left_Up, 'replicate', 'same', 'corr');
  I_RU = imfilter(I, K_Right_Up, 'replicate', 'same', 'corr');
  I_C = imfilter(I, K_Center, 'replicate', 'same', 'corr');
  I_LD = imfilter(I, K_Left_Down, 'replicate', 'same', 'corr');
  I_RD = imfilter(I, K_Right_Down, 'replicate', 'same', 'corr');
  M_LU = (I>(I_LU+mu));
  M_RU = (I>(I_RU+mu));
  M_C = (I>(I_C+mu));
  M_LD = (I>(I_LD+mu));
  M_RD = (I>(I_RD+mu));
  M = M_LU.*M_RU.*M_C.*M_LD.*M_RD;
  Mark = M(:, :, 1).*M(:, :, 2).*M(:, :, 3);
  Mark = 1-Mark;
  % U-V space
  C = 1/3.*(I(:, :, 1) + I(:, :, 2) + I(:, :, 3));
  U = (2*C-I(:, :, 2)-I(:, :, 3))./C;
  V = max(1-I(:, :, 2)./C, 1-I(:, :, 3)./C);
  dis = U.*U + V.*V;
  M = (dis>epsilon);
  Mark = max(M, Mark);
end