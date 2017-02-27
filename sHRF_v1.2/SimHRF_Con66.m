function [h] = SimHRF_Con66(V)
% 
% function [h] = SimHRF_Con55(V)
% 
% This function generate HRF using 2 gamma function ( 6 paramters) with
% legth of 30 seconds. 
% V = parameters
% h = HRF time course
%
% Date: Sep 23, 2011
%
% Author: Zuyao Shan

j = 0;
for t = (0:0.1:29.9)
    j=j+1;
    h(j,2)= V(1) * (t^(V(2)-1)) * (V(3)^V(2)) * (exp(-1*V(3)*t))/gamma(V(2)) ...
                    - V(4) * (t^(V(5)-1)) * (V(6)^V(5)) * (exp(-1*V(6)*t))/gamma(V(5));
    h(j,1)= t;
end
end

