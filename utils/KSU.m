function [rst] = KSU(vector,penalization)
%KSU Summary of this function goes here
%   Detailed explanation goes here
N = length(vector);
rst = 1/penalization * (1/N*sum(exp(penalization*vector)));
end

