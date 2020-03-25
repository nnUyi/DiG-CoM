function [SE, w1, w2] = RotateFilter(theta)
Grax = [
    0 0 0 
    0 1 0
    -1 0 0];
% Gray2 = [0 0 0;0 1 0;0 0 -1];
Gray = [
    0 0 0
    0 1 0
    0 -1 0];
tmp1 = 1/sin(pi/4-theta);
tmp2 = 1/sin(theta);
% w1 = exp(pi/8-theta);
% w2 = exp(theta-pi/8);
w2 = tmp1/(tmp1+tmp2);
w1 = tmp2/(tmp1+tmp2);
SE = w1.*Grax+w2.*Gray;
SE = SE./(w1+w2);
% SE = (pi/4-theta).*Grax+theta.*Gray;
% SE = SE./(pi/4);
end