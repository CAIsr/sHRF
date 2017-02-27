function [ Sti_tc ] = sSti_Stick(onsets,T,m)
%
%   Function [Sti_tc] = sSti_Stick(onsets,T)
%
%   Generate a time series of box car function 
%       Sti_tc: n by 2 matrix with first column of time and second column of
%                values
%       onsets: Times of onset 1 by n array
%       T: totoal time length of the session
%       
%   Author: Zuyao Shan
%   
%   Date: Sep 23, 2011
%
    Times = (0:0.1:(T-0.1));
    Sti = zeros(1, size(Times,2));
    for i=(1:1:size(onsets,2))
        Sti(round(onsets(i)*10)) =1/m;        
    end
    Sti_tc=[Times; Sti];
    Sti_tc =Sti_tc';    
end

