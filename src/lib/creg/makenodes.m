function nodeStruct = makenodes(p,t,order)

%MAKENODES produces:
%  a matrix NODES containing coordinates for all of thenodes to be used, 
%  a matrix NODEINDEX defining which nodes correspond to each element.  
%  If NORDER is 2, the midpoint of each edge is computed and added to 
%  POINTS to obtain matrix NODES. 
%  The row index of that midpoint is then added to the rows of TRIANGLES 
%  containing that edge to define NODEINDEX.
%  If NORDER is 1, nodes corresponds to vertices and NODEINDEX is 
%  identical to TRIANGLES.
%EXT is a vector containing the coordinates of the exterior nodes of 
%the mesh.
%
%  nvert:  number of vertices
%  nele:   number of triangles or elements
%
%Input: POINTS is an nvert by 2 matrix containing the x and y
%	coordinates of each of the nvert points in the right-hand rule mesh.
%   POINTS = P' where P is the points matrix for pde
%   The call can use P directly (see below).
%	TRIANGLES is T(1:3,:)' where T is the triangle index matrix from pde.
%   Vertices must be numbered in the counterclockwise direction.
%   NORDER is the order of elements, and may be either 1 or 2 (default)
%
%Output: NODES:  a nnodes*2 matrix whose i'th row contains
%		the coordinates of the i'th nodal variable.
%       Nodes for the second order element consist of vertices and 
%       midpoints of edges, that is, 6 per triangle.
%		The first NVER rows of NODES is POINTS, and the remainder are
%       the edge midpoints.
%       Nodes for the first order element consist of only vertices.
%
%
%	 NODEINDEX:  for NORDER == 2, an nele*6 matrix whose i'th row
%		contains the row numbers (in NODES) of the
%		nodal variables defining the i'th finite 
%		element.  If the i'th row of FMESH is [V1 V2 V3]
%		then the i'th row of nodeindex is
%		[V1 V(12) V2 V(23) V3 V(31)], where Vi is the
%		row number of the i'th point and V(ij) is the 
%		row number of the midpoint of the edge defined
%		by the i'th and j'th points.
%       If NORDER == 1, NODEINDEX is TRIANGLES.
%
% Last modified on August 19, 2010 by Jim Ramsay.

%  Set default order

if nargin < 3, order = 2;  end

if size(p,2)>2
    %  transpose if points and triangles come from pde toolbox
    nodes = p';
    t = t(1:3,:)';
else
    nodes = p;
end;

nele = size(t,1);  %  number of elements
nver = size(p,1);     %  number of vertices

Jvec   = zeros(nele,1);    %  vector of jacobian values
metric = zeros(nele,2,2);  %  3-d array of metric matrices

if order == 2
    
    rec = sparse(nver,nver);
    ind = [ 1 2 ; 2 3 ; 3 1 ];
    nodeindex = zeros(nele,6);
    nodeindex(:,[1 3 5]) = t;
    
    for i = 1:nele
        for j = 1:3
            if rec(t(i,ind(j,1)),t(i,ind(j,2)))==0
                nodes = [nodes; [.5 .5]*nodes(t(i,ind(j,:)),:)];
                rec(t(i,ind(j,1)),t(i,ind(j,2))) = ...
                    size(nodes,1);
                rec(t(i,ind(j,2)),t(i,ind(j,1))) = ...
                    size(nodes,1);
                nodeindex(i,2*j) = size(nodes,1);
            else
                nodeindex(i,2*j) = ...
                    rec(t(i,ind(j,1)),t(i,ind(j,2)));
            end;
        end;
        
        %  deviations of vertices 2 and 3 from vertex 1
        
        diff1x = nodes(nodeindex(i,3),1) - nodes(nodeindex(i,1),1);
        diff1y = nodes(nodeindex(i,3),2) - nodes(nodeindex(i,1),2);
        diff2x = nodes(nodeindex(i,5),1) - nodes(nodeindex(i,1),1);
        diff2y = nodes(nodeindex(i,5),2) - nodes(nodeindex(i,1),2);
        
        %  Jacobian or area of triangle
        
        Jvec(i) = diff1x.*diff2y - diff2x.*diff1y;
        
        %  Compute contravariant transformation matrix
        
        Ael = [ diff2y -diff1y; ...
               -diff2x  diff1x]./Jvec(i);

        %  Compute metric matrix
        
        metric(i,:,:) = Ael'*Ael;
    end;
elseif order == 1
    nodeindex = t(:,1:3);
    
    for i=1:nele
        
        %  deviations of vertices 2 and 3 from vertex 1
        
        diff1x = nodes(nodeindex(i,2),1)-nodes(nodeindex(i,1),1);
        diff1y = nodes(nodeindex(i,2),2)-nodes(nodeindex(i,1),2);
        diff2y = nodes(nodeindex(i,3),2)-nodes(nodeindex(i,1),2);
        diff2x = nodes(nodeindex(i,3),1)-nodes(nodeindex(i,1),1);
        %  Jacobian or area of triangle
        
        Jvec(i) = diff1x.*diff2y - diff2x.*diff1y;
        
        %  Compute contravariant transformation matrix
        
        Ael = [ diff2y -diff1y; ...
               -diff2x  diff1x]./Jvec(i);

        %  Compute metric matrix
        
        metric(i,:,:) = Ael'*Ael;
    end
else
    error('ORDER not 1 or 2.');
end

nodeStruct.order     = order;
nodeStruct.nodes     = nodes;
nodeStruct.nodeindex = nodeindex;
nodeStruct.J         = Jvec;
nodeStruct.metric    = metric;