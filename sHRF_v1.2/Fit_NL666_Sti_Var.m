function [hrf, fit, e, param, aic, VM] = Fit_NL666_Sti_Var(TC,STI)
% function [hrf, fit, e, param] = Fit_NL(TC,STI)
%
% Fit the time course with three summed cannocial functions using simplex 
%   method  
%
% INPUTS:
% 
% TC    - time course
% STI   - time course of stimulus function
% T     - estimated HRF function length
%  
% OUTPUTS:
%
% hrf   - estimated hemodynamic response function
% fit   - estimated time course
% e     - residual errors
% param - estimated amplitude, height and width
%
% Author:Zuyao Shan 
% 
% Date: Aug 10, 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit the two gamma model

%%%%%%%%%% Default parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0 = [6 7 1 1 16 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];        % initial values for gamma functions
n_V= 6;     % number of parameters


LB = [0 2 0.5 0 6 0 -20 0 0 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9 -0.9];  %set the search lower bound 
UB = [10 12 2 6 25 1.5 20 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9]; %set the search upper bound

% Find optimal values
options = optimset('MaxFunEvals',10000,'Maxiter',10000,'TolX',1e-10,'TolFun',1e-10,'Display','off');
VM = fminsearchbnd(@msq_gamma,V0,LB,UB,options,TC,STI);

% Use optimal values to fit hemodynamic response functions
hrf = SimHRF_Con66(VM(1, 1:6));                   % Calculate HRF estimate (fit, given theta)

param = zeros(5);
fit = zeros(size(TC,1),2);
e =0;
aic=0;

param = HRF_para(hrf);
fit = conv((STI(:,2))',(hrf(:,2))');
fit = fit + VM(1,7);
for i=(1:1:size(TC,1))
    e = e +(TC(i,2)-fit(floor(TC(i,1)*10+1))) * (TC(i,2)-fit(floor(TC(i,1)*10+1)));
end
aic = size(TC,1).*log(e/size(TC,1))+ (2 .* n_V .* (n_V + 1))./(size(TC,1)- n_V -1);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SUBFUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m=msq_gamma(V,TC,STI)
h = SimHRF_Con66(V(1,1:6));% Get gamma model corresponding to parameters V
STI = Sti_Var(STI,V(1,8:24));
yhat = conv(STI(:,2)',h(:,2)');       % Convolve gamma model with stimulus function
% Calculate cost function
m = 0;
for i=(1:1:size(TC,1))
    m = m +abs(TC(i,2)-yhat(floor(TC(i,1)*10+1))-V(1,7)); %* (TC(i,2)-yhat(floor(TC(i,1)*10+1)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


