function [Vertices Faces FaceID] = mercator(n, shapetype)
% function [Vertices Faces FaceID] = MERCATOR(n, shapetype))
%
% Generate a cartesian coordinates of the sphere of radius 1 mapped
% from Mercator (cylindrical longitude/latitude)
%
% INPUT:
%   n: 2 x 1 array: number of linear subdivision in latitude/longitude
%   shapetype: optional 'spherical' (default) or 'spheroid' (earth shape)
% OUTPUT:
%   Vertices: (nv x 3) array, where nv = (n(1)+1)*n(2), 3D vertices coordinates
%   Faces:    (nf x 4) array, where nf = n(1)*n(2), vertices indices of patches
%   FaceID:   (nf x 1) array, number where the faces belong (to the cube
%                             topology, 1)
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%   original: 19-Aug-2009
%   20-Aug-2009: spheroid shape

n = round(n(:).');
if any(n)<1
    error('MERCATOR: n must be stricly positive number');
end

classIdx = 'uint32';
if double(intmax(classIdx)) < prod(n+1)
    classIdx = 'double'; % You will have anyway memory problem later
end

if nargin<2 || isempty(shapetype)
    shapetype = 'spheroid';
end

%% Discretization and mercator projection
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

lat = linspace(-pi/2, pi/2, n(1)+1);
lon = linspace(-pi, pi, n(2)+1);
lon(end)=[];
[LON LAT] = meshgrid(lon, lat);
X = cos(LON).*cos(LAT);
Y = sin(LON).*cos(LAT);
Z = rr*sin(LAT);

%% Generate face
Vertices = [X(:) Y(:) Z(:)];
clear X Y Z

%% Patch topology of the first face
vx = feval(classIdx, 1:n(2));
vy = feval(classIdx, 1:n(1));
[i j] = ndgrid(vy,vx);
i = i(:); j = j(:);
vxp1 = vx+1; vxp1(end) = 1;
vyp1 = vy+1;
[ip1 jp1] = ndgrid(vyp1,vxp1);
ip1 = ip1(:); jp1 = jp1(:);
% i-> vertical, j-> horizontal
Faces = [sub2ind([n(1)+1 n(2)], i, j) ...
         sub2ind([n(1)+1 n(2)], i, jp1) ...
         sub2ind([n(1)+1 n(2)], ip1, jp1) ...
         sub2ind([n(1)+1 n(2)], ip1, j)];

% Faces ID
FaceID = repmat(uint8(1),prod(n),1);
FaceID = FaceID(:);

end % MERCATOR

