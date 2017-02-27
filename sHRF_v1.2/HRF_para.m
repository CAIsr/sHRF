function [ param ] = HRF_para( hrf )
%
%function [param] = HRF_para(hrf)
%Estimate the HRF parameters of height, time to peak, and width of 
%positive response  
%   param: estimated HRF parameters of positive response
%   hrf: n by 2 hrf matrix with first column of time and seconf column of
%       HRF values
%   
%Author: Zuyao Shan
%
%Date: August 9, 2011
[MaxTab(1),MaxTab(2)]= max(hrf(:,2));
param(1)= MaxTab(1);
param(2)= hrf(MaxTab(2),1);
Half_Heights = intercept((hrf(:,2))', (MaxTab(1) ./ 2));
if size(Half_Heights,1)>1
    param(3) = hrf(Half_Heights(2),1)-hrf(Half_Heights(1),1);
else
    param(3) = NaN;
end
param(4)=0;
if size(Half_Heights,1)>1
    for i=(1:1:Half_Heights(1))
        if ((hrf(i,2)-0.1*param(1)) > 0) && param(4)==0
            param(4)=hrf(i,1);
        end
    end
end
%%
%*****Calculate the area under the curve***********************************
lc = MaxTab(2);
a_lc = 0;
while ((lc > 1) && ((hrf(lc,2)-0.0001) > 0))
    a_lc = a_lc + hrf(lc,2) * 0.1;
    lc = lc - 1;
end
rc = MaxTab(2)+1;
a_rc = 0;
while ((rc < 300) && ((hrf(rc,2)-0.0001) > 0))
    a_rc = a_rc + hrf(rc,2) * 0.1;
    rc = rc + 1;
end
param(5)=a_lc+a_rc;
%%
%Parameter check to exclude non sensible value, if outboud data put as NaN
if param(1) == 0
    param(1:5) = NaN;
end
if param(2) > 10
    param(1:4) = NaN;
end
if param(1) < 0.05
    if param(2) > 8
       param(1:5) = NaN;
    end
    if param(3) > 10
       param(1:5) = NaN;
    end
end
if param(2) < 1
    if param(1) > 4
      param(1:5) = NaN;
    end
end
end

