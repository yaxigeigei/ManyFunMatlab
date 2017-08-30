function [  ] = AddXunit(unit)
%ADDXUNIT Summary of this function goes here
%   Detailed explanation goes here


xlb = get(gca, 'XTickLabel');
l = size(xlb, 2);
xlb(end, l+1:l+length(unit)) =  unit;
set(gca, 'XTickLabel', xlb);



end

