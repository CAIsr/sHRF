function [ ValueTab ] = intercept(V, x)
%
% function [ValueTab] = intercept(V, x)
%
% Find the elements in array at which its value equal to x
%   ValueTab: an array of positions where the array have the
%             same/approximate value of x
%   V: input array
%   x: scalar value 
%   
% Author: Zuyao Shan
%
% Date: Aug 5, 2011
ValueTab=[];

Temp_V = V - x;
size_V = size(V);

for i = (1:1:size_V(2)-1)
    j = i+1;
    if Temp_V(i) > 0 && Temp_V(j) < 0
        ValueTab = [ValueTab; i];
    end
    if Temp_V(i) < 0 && Temp_V(j) > 0
        ValueTab = [ValueTab; i];
    end
end
end

