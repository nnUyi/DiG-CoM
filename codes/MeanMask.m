function[F] = MeanMask(Input, w1, w2, Mask)
%%       Mean filtering with a introduced Mask
% to improve the afficiency

filter = ones(w1, w2);
I_temp = Input.*Mask;
F = imfilter(I_temp, filter, 'replicate', 'same');
filM = imfilter(Mask, filter, 'replicate', 'same');
F = F./filM;

tmpf = F;
tmpf(isnan(tmpf)) = Input(isnan(tmpf));
F = tmpf;

tmpf = F;
tmpf(tmpf==Inf) = Input(tmpf==Inf);
F = tmpf;
end