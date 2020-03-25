function [F] = BilateralMark(Input, w, Mark, sigma_s, sigma_r, T)
%%                       Parameters
  % F: image after bilateral filters
  % Input: input image, 3 channels: R G B, double range 0 and 1
  % w: window size of filters, should be odd integer
  % sigma_s: gaussian parameter for spatial filters
  % sigma_r: gaussian parameter for range filters
  % T: 255 or 1
% Input = round(255*Input);
gamma = pi/(2*T);
rho = gamma*sigma_r;
if rho >= 1
    N = 15;
else
    N = ceil(1/(rho*rho));
end
[m, n, c] = size(Input);
F = zeros(m, n, c); 
tmp = rho*sqrt(N);
for channel=1:c
    B = zeros(m, n);  % B: sum of bn(x)
    G = zeros(m, n);  % G: sum of gn(x)
    for num=0:N
        I = Input(:, :, channel);
        I_temp = I.*Mark;
        index = gamma*(2*num-N)/tmp;
        bn = exp(1i*index*I).*Mark;
        gn = exp(1i*index*I_temp).*I_temp;
        dn = 2^(-N)*nchoosek(N, num)*exp(-1i*index*I);
        %%            Filter with gaussian spatial kernel
        step = (w-1)/2;
        [X, Y] = meshgrid(-step:step, -step:step);
        Gau = exp(-(X.^2+Y.^2)/(2*sigma_s^2));
        % Gau = 1/sum(Gau(:)).*Gau;
        bn_par = imfilter(bn, Gau, 'same', 'symmetric', 'conv');
        gn_par = imfilter(gn, Gau, 'same', 'symmetric', 'conv');
        B = B + bn_par.*dn;
        G = G + gn_par.*dn;
    end
    F(:, :, channel) = real(G./B);
   tmpf = F(:, :, channel);
   tmpf(isnan(tmpf)==1) = I(isnan(tmpf)==1);
   F(:, :, channel) = tmpf;
end
% F = double(F/255);
end