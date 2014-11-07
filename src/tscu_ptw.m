function [distance,path1,path2] = tscu_ptw(x,y,~)
%TSCU_PTW Parametric Time Warping
%   [DISTANCE PATH1 PATH2]=TSCU_PTW(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
[w,~] = quadwarp(x', y');
%fprintf('w: %f\n',w);
%fprintf('sel: %f\n',sel);
w(w<1)=1;
w(w>length(x)) = length(x);

path1=w';
path2=1:length(y);


xAligned = interp1(1:length(x),x,path1);
yAligned = interp1(1:length(y),y,path2);
distance=norm(xAligned-yAligned);
if isnan(distance)
   fprintf('%10.5f %10.5f\n',w(1),w(end));
end
end