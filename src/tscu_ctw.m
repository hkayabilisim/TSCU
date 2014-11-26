function [distance,path1,path2] = tscu_ctw(x,y,~)
%TSCU_CTW Canonical Time Warping
%   [DISTANCE PATH1 PATH2]=TSCU_CTW(X,Y,OPTIONS) aligns
%   the objects X and Y using available options OPTIONS and returns
%   the distance in DISTANCE and warping paths in PATH1 and PATH2.
aliUtw = utw({x,y}, [], []);
aliCtw = ctw({x,y}, aliUtw, [], [], [], []);

path1 = aliCtw.P(:,1)';
path2 = aliCtw.P(:,2)';

xAligned = interp1(1:length(x),x,path1);
yAligned = interp1(1:length(y),y,path2);
distance=norm(xAligned-yAligned);
end