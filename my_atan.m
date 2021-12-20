function value = my_atan(x1,y1,x2,y2)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
if x1<x2 && y1<y2 || x1<x2 && y1>=y2
    value = atan((y2-y1)/(x2-x1));
elseif x1>=x2 && y1<y2
    value = atan((y2-y1)/(x2-x1))+pi;
elseif x1>=x2 && y1>=y2
    value = atan((y2-y1)/(x2-x1))-pi;
end
end