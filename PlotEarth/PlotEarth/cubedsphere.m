function [Vertices, Faces, FaceID] = cubedsphere(n, prjtype, shapetype)
% function [Vertices Faces FaceID] = cubedsphere(n, prjtype, shapetype))
%
% Generate a cartesian coordinates of the cubed-sphere of radius 1
% the grid is gnomonic equiangular/equidistance central projection
%
% INPUT:
%   n: number of linear subdivision of the face
%   prjtype: optional, 'equiangular' (default) or 'equidistance'
%   shapetype: optional 'spherical' (default) or 'spheroid' (earth shape)
% OUTPUT:
%   Vertices: (nv x 3) array, where nv = n^2*6+2, 3D vertices coordinates
%   Faces:    (nf x 4) array, where nf = n^2*6, vertices indices of patches
%   FaceID:   (nf x 1) array, number where the faces belong (to the cube
%                             topology, 1-6)
%      +---+
%      | 6 |
%  +---+---+---+---+
%  | 4 | 1 | 2 | 3 |
%  +---+---+---+---+
%      | 5 |
%      +---+
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%   14-Aug-2009: original
%   15-Aug-2009: change the patch order (ndgrid -> meshgrid, code #131)
%   17-Aug-2009: return Faces UINT32 and FaceID as UINT8 class
%                invoke UNIQUE for boundary pixels only
%   18-Aug-2009: cosmetic changes
%   20-Aug-2009: spheroid shape
%   12-Jul-2014: Fix casting pb in 2014

n = round(n);
if n<1
    error('CUBEDSPHERE: n must be stricly positive number');
end

if nargin<2 || isempty(prjtype)
    prjtype = 'equiangular';
end

if nargin<3 || isempty(shapetype)
    shapetype = 'spheroid';
end

n = n+1;

classIdx = 'uint32';
if double(intmax(classIdx)) < n*n*6
    classIdx = 'double'; % You will have anyway memory problem later
end

%% Discretize basic face of a cube
switch lower(prjtype)
    case 'equiangular'
        % equidistance projection
        x = 1;
        theta = linspace(-pi/4, pi/4, n);
        y = tan(theta);
        z = y;
        [X, Y, Z] = ndgrid(x,y,z);
    case 'equidistance'
        % equidistance projection
        x = 1;
        y = linspace(-1,1,n);
        z = y;
        [X, Y, Z] = ndgrid(x,y,z);
    otherwise
        error('CUBEDSPHERE: prjtype must be ''equiangular'' or ''equidistance''');
end

%% Shape
switch lower(shapetype)
    case 'spheroid'
        % Radius of the earth; equatorial/polar
        % http://en.wikipedia.org/wiki/Earth_radius
        requator = 6378.1370;
        rpolar = 6356.7523;
        rr = rpolar/requator;
    case 'spherical'
        rr = 1;
    otherwise
        error('Mercator: unknown shapetype');
end

%% Project on sphere S2
S = sqrt(1./(X.^2+Y.^2+Z.^2));
X = X.*S;
Y = Y.*S;
Z = rr*Z.*S;

%% Generate six faces
C = zeros([numel(X) 6 3],'double');
C1 = [X(:) Y(:) Z(:)];
MG = makehgtform('zrotate',0);
MG = MG(1:3,1:3);
% Generate six faces
C(:,1,:) = C1*(MG.');
% Other faces are rotations of the first
M = makehgtform('zrotate',pi/2);
C(:,2,:) = C1*((MG*M(1:3,1:3)).');
M = makehgtform('zrotate',pi);
C(:,3,:) = C1*((MG*M(1:3,1:3)).');
M = makehgtform('zrotate',-pi/2);
C(:,4,:) = C1*((MG*M(1:3,1:3)).');
M = makehgtform('yrotate',pi/2);
C(:,5,:) = C1*((MG*M(1:3,1:3)).');
M = makehgtform('yrotate',-pi/2);
C(:,6,:) = C1*((MG*M(1:3,1:3)).');
C = reshape(C, [], 3);

clear X Y Z C1 S

%% Detect common vertices (between two neighboring faces)
% Index of points on the boundary 
ileft = feval(classIdx,1:n);
itop = feval(classIdx,1:n:n^2);
ibottom = feval(classIdx,n:n:n^2);
iright = feval(classIdx,(n-1)*n+1:n^2);
iboundary = munion(ileft, iright, itop, ibottom);
offset = feval(classIdx, (n^2)*(0:5));
%iboundary = bsxfun(@plus, offset, iboundary(:));
iboundary = repmat(offset, length(iboundary), 1) + ...
            repmat(iboundary(:), 1, length(offset));
iboundary = sort(iboundary(:));

tol = 1e-3/(n-1); % tolerance below which two pixels are considerd as equal
Cbr = round(C(iboundary,:)/tol);
[trash Is Js] = unique(Cbr,'rows'); %#ok mlint complains on trash

m = size(C,1);
duplicate = true(size(iboundary));
duplicate(Is) = false;

% Remove index of duplicate vertices
I = feval(classIdx,(1:m).');
I(iboundary(duplicate)) = [];

% Reverse index mapping: J(k) is the new index from the old index k
J = ones(m,1,'double');
J(iboundary(duplicate)) = 0;
J = feval(classIdx,cumsum(J));
% populate the iboundary index
selected = J(iboundary(Is));
J(iboundary) = selected(Js);

% clean up
clear duplicate Cbr trash selected Is Js

% Remove vertices
Vertices = C(I,:);

clear C I
%% Patch topology of the first face
v = feval(classIdx, 1:n-1);
[i, j] = meshgrid(v,v); % 090815: Change to meshgrid from ndgrid
i = i(:); j=j(:);
% i-> horizontal, j-> vertical
Faces = [sub2ind([n n], i, j) ...
         sub2ind([n n], i+1, j) ...
         sub2ind([n n], i+1, j+1) ...
         sub2ind([n n], i, j+1)]; % sub2ind cast the result to double since V2014
clear i j 
% Topology for all faces
Faces = reshape(Faces,[size(Faces,1) 1 size(Faces,2)]);
Faces = feval(classIdx, Faces); % sub2ind cast the result to double since V2014
offset = feval(classIdx, (n^2)*(0:5));
Faces = bsxfun(@plus, offset, Faces);
% Faces = repmat(offset, [size(Faces,1) 1 size(Faces,3)]) + ...
%         repmat(Faces, [1 length(offset) 1]);
Faces = reshape(Faces, [], 4);

% update the face topology
for k=1:size(Faces,2) % use for-loop to avoid create temporary large array
    Faces(:,k) = J(Faces(:,k));
end
clear J

% Faces ID
FaceID = uint8(1:6);
FaceID = repmat(FaceID,(n-1)^2,1);
FaceID = FaceID(:);

end % CUBEDSPHERE

function u = munion(varargin)
% union of a list of index
u = [];
for k=1:length(varargin)
    u = union(u, varargin{k});
end
end % munion
