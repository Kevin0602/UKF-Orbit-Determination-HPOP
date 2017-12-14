function h = plotearth(varargin)
% function h = plotearth(Property1, Value1, ...)
% Plot the earth map on sphere
%
% Propert/Value pairs, Properties are strings:
% 'MapType', string map among the following
%   - low resolution: 'texture' (default), 'relief'
%   - hight resolution: 'nasa', 'avhrr', 'usgs' (large memory required)
%   - NEO maps, hight resolution: 'bluemarble' 'topo'
% 'NEOMap': filename of the maps from NEO/ESA website (see below)
% 'SampleStep', integer: 1 -> highest resolution from the orginal map
%                      larger -> lower resolution
% 'FullColor', boolean: true -> RGB (required more memory)
%                       false -> indexed color (degraded color, default)
% 'NColors', integer: number of color used for indexed color map
%                     256 by default. Ignored when 'FullColor' is TRUE
% 'ViewPoint', 1 x 2 array, azimuthal/elevation angles (degree)
% 'Shape': 'spheroid' (default), or 'spherical'
% 'Plot2DMap', boolean: plot the resampling flat 2D map, false by default
% 'AxeHandle', handle of the axes
%
% Some plots is based on the maps available from
%   http://www.progonos.com/furuti/MapProj/CartIndex/cartIndex.html
%
% Few other are based from Nasa Earth Observation (NEO) 
%   http://neo.sci.gsfc.nasa.gov/Search.html
% Instruction to download NEO map from this website:
%   - Select your map
%   - Select option <Full>, <Color>, <PNG> on the right panel
%   - Click on <GET IMAGE> button and save (right mouse click) 'mymap.png'
%     on the same directory, filename 'Population.png'
%   - Under Matlab
%       >> plotearth('neomap', 'Population', ...)
%
% See the respective websites for copyright information
%
% A collection of maps is compiled and can be downoaded in the link
%   http://www.mediafire.com/?m2mn1mgdngt
%
% OUTPUT: Return the handle of the AXE
%
% NOTES:
%   - Shape is spheroid by default, unless id SPHERICAL is specified 
%   - The earth center is located at (0,0,0)
%   - Greenwitch meridian is aligned with the plane y=0
%   - Units are normalized by the eath equatorial radius of 6378.1370 km
%
% See also: cubedsphere, mercator
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% History:
%   17-Aug-2009, original
%   18-Aug-2009, better use of memory
%                correct BUG for color-indexing
%                add 'NColors' property
%                indexing color is now default option
%                small gaps between cube-faces are now removed
%   19-Aug-2009, NEO/ESA map
%                correct bug overflow colormap in original file
%   20-Aug-2009, Read tiff/jpg maps
%                Spheroid shape
%   18-Oct-2011, fix Plot2DMap lowercase bug 

%% Parse inputs

% Default values
inputs = struct('maptype', 'texture', ...
                'neomap', '', ...
                'samplestep', NaN, ...
                'fullcolor', false, ...
                'verbose', true, ...
                'plot2dmap', false, ...
                'axehandle', NaN, ...
                'ncolors', 256, ...
                'viewpoint', [155 15],...
                'shape','spheroid');
            
validprop = fieldnames(inputs);
if mod(length(varargin),2)
    error('PLOTEARTH must be called with property/value pair');
end
% Parsing inputs
for k=1:2:length(varargin)
    property = varargin{k};
    if ~ischar(property)
        error('PLOTEARTH property must be a string');
    end
    property = lower(strtrim(property));
    idx = strmatch(property,validprop);
    if isempty(idx)
        error('PLOTEARTH not valid property ''%s''', property);
    end
    property = validprop{idx(1)};
    value = varargin{k+1};
    if ~isempty(value)
        inputs.(property) = value;
    end
end

% Displatch the inputs
maptype = inputs.maptype;
neomap = inputs.neomap;
samplestep = round(inputs.samplestep);
fullcolor = inputs.fullcolor;
verbose = inputs.verbose;
plot2dmap = inputs.plot2dmap;
axehandle = inputs.axehandle;
ncolors = inputs.ncolors;
viewpoint = inputs.viewpoint;
shape = inputs.shape;

%% Map selection
if ~isempty(neomap) % NEO map
    mapfname = neomap;
    % Add the png extension
    [pname maptype ext] = fileparts(mapfname); %#ok
    if isempty(ext)
        % Automatically look for image with following extension
        for ext = {'.png' '.PNG' ...
                   '.tif' '.TIF' '.tiff' '.TIFF' ...
                   '.jpg' '.JPG' '.jpeg' '.JPEG'}
            if exist([mapfname ext{1}],'file')
                mapfname = [mapfname ext{1}]; %#ok
                break;
            end
        end
    end
    projection = 'mercator';
    offsetx = 1; offsety = 1;
    nxy = NaN;
    mstep = 1;
else
    switch lower(maptype)
        case 'relief'
            mapfname = 'cbGn-s100.rel.png';
            projection = 'cube';
            offsetx = 2; offsety = 27; % upper/left corner
            nxy = 200;                 % size of each cube face, in pixel
            mstep = 1;                 % sampling step
        case 'texture'
            mapfname = 'cbGn-s100.tex.png';
            projection = 'cube';
            offsetx = 3; offsety = 27;
            nxy = 200;
            mstep = 1;
        case 'nasa'
            mapfname = 'cbGn_pof-bm.png';
            projection = 'cube';
            offsetx = 465; offsety = 215;
            nxy = 1524;
            mstep = 2;
        case 'avhrr'
            mapfname = 'cbGn_pof-z-25-pat.png';
            projection = 'cube';
            offsetx = 465; offsety = 215;
            nxy = 1524;
            mstep = 2;
        case 'usgs'
            mapfname = 'cbGn_poe-z20-tex.png';
            projection = 'cube';
            offsetx = 466; offsety = 210;
            nxy = 1525;
            mstep = 2;
        case 'bluemarble'
            mapfname = 'BlueMarble.png';
            projection = 'mercator';
            offsetx = 1; offsety = 1;
            nxy = NaN;
            mstep = 1;
        case 'topo'
            mapfname = 'Topo.png';
            projection = 'mercator';
            offsetx = 1; offsety = 1;
            nxy = NaN;
            mstep = 1;
        otherwise
            fprintf(['Maptype should be: ' ...
                '''texture'' ''relief'' ''nasa'' ''avhrr'' ''usgs''\n']);
            error('PLOTEARTH: unknown maptype <%s>', maptype);
    end
end

% Use default samplestep associated to each map 
if isnan(samplestep)
    samplestep = mstep;
end

    % nested function to display a text on command line window
    function display(varargin)
        if verbose
            fprintf(varargin{:});
        end
    end

%% Read the map
display('read the map\n')
display('\tmap: %s\n', mapfname)
[A map] = imread(mapfname);

try % sometime memory problem occurs
    
    % rotate 90 degree
    if strcmp(projection,'cube')
        A = permute(A,[2 1 3]);
        A = flipdim(A,1);
    end
    
    % Map to RGB
    if ~isempty(map)
        szA = size(A);
        if ~strcmp(class(A),'double') % must do that to avoid overflow!
            A = single(A); % this should cover 1/2-byte classes
        end
        A=map(A+1,:);
        A=reshape(A,[szA 3]);
        clear map
    end

    % Resampling RGB
    display('resampling image\n');
    [A m n] = resampling(A, offsetx, offsety, nxy, samplestep, ...
                         projection, @display);
    display('\tresolution = %d x %d x %d\n', m, n(1), n(2));
  
    % Indexed color
    if ~fullcolor && size(A,3)==3
        display('indexing color\n');
        if isempty(which('rgb2ind'))
            warning('RGB2IND:Missing', ...
                    'RRB2IND is missing (Image processing required)');
        else
            [A clmap] = rgb2ind(A,ncolors);
            display('\tnumber of colors = %d\n', size(clmap,1));
        end
    else
        display('full RGB color scheme\n');
    end
    % 3 (RGB) or 1 (indexed color)
    d3 = size(A,3);
    
    % Plot the flat 2D map
    if plot2dmap
        fig = figure();
        set(fig, 'Name', ['Map: ' maptype]);
        ax=axes('Parent',fig);
        image(A,'Parent',ax);
        if exist('clmap','var')
            colormap(ax,clmap);
        end
        axis(ax,'equal');
    end
    
    % Project on sphere
    switch projection
        case 'cube',
            display('compute cubed-sphere\n');
            [Vertices Faces FaceID] = cubedsphere(n(1), ...
                                               'equidistance', shape);
        case 'mercator',
            display('compute mercator-sphere\n');
            [Vertices Faces FaceID] = mercator(n, shape);            
    end
    
    %% Preparation
    display('preparation\n')
    
    xleft = 1 + ((1:m)-1)*n(2);
    xright = xleft + (n(2)-1);
    yupper = 1 + ((1:1)-1)*n(1);
    ylower = yupper + (n(1)-1);
    
    % Reshape the map as long vector of RGB
    xl = xleft(1:m);
    xr = xright(1:m);
    yu = yupper(repmat(1,m,1));
    yl = ylower(repmat(1,m,1));   
    
    % Save half of memory for full-color plot
    if d3==3
        A = single(A);
    end
    
    %% Plot the earth map on sphere
    display('plot earth\n')
    % Open axe
    if ~ishandle(axehandle)
        fig = figure();
        set(fig, 'Color','k', 'Name', ['Earth: ' maptype]); % black background
        
        pos = get(0,'screensize'); % full screen
        set(fig, 'Position', pos);
        ax = axes('Parent',fig,'pos',[0 0 1 1]);
    else
        ax = axehandle;
    end
    hold(ax,'on');
    
    % Mapping method
    if d3==3
        Mapping = 'direct';
    else
        Mapping = 'scaled';
    end
    
    for ID=1:m % loop on faces   
        
        % Get colors
        x = xl(ID):+1:xr(ID);
        y = yl(ID):-1:yu(ID);
        subs = {y x ':'};
        CData = reshape(A(subs{:}),[n(1),n(2),d3]);
        if ~strcmp(CData,'uint8')
            CData = double(CData); % CData must be double or uint8
        end
        
        % Get grid coordinates from vertices
        IDX = GridIdx(Faces, FaceID, ID, n, projection);
        x = reshape(Vertices(IDX,1),[n(1)+1 n(2)+1]);
        y = reshape(Vertices(IDX,2),[n(1)+1 n(2)+1]);
        z = reshape(Vertices(IDX,3),[n(1)+1 n(2)+1]);
        
        % Plot
        surf(ax,x,y,z,...
               'CData',CData,... 
               'EdgeColor','none',...
               'CDataMapping',Mapping);
    end % for loop
    if exist('clmap','var')
        colormap(ax,clmap);
    end    
    axis(ax,'off');
    axis(ax,'equal');
    view(ax,viewpoint);
    
    if nargout>=1
        % return the handle
        h = ax;
    end
    
catch %#ok % 'catch ME': syntax not compatible with V2006
    ME = lasterror(); %#ok
    if strcmp(ME.identifier,'MATLAB:nomem')
        fprintf('PLOTEARTH: Memory is full\n');
        fprintf('\tClose you figures (if any) and clear workspace\n');
        fprintf('\tIf it''s still fail, use a large SAMPLESTEP (=%d here)\n', samplestep);
    end
    rethrow(ME);
end

end % plotearth

%%
function [Acoarse m n] = resampling(A, x1, y1, nxy, samplestep, ...
                                    projection, displayfun)
% function [Acoarse m n] = resampling(A, x1, y1, nxy, samplestep, ...
%                                    projection, displayfun)
% Resampling the original in coarser resolution
% INPUTS:
%   (x1,y1) is the index of the upper/left corner
%   nxy is the linear size of the rectangular face (pixels)
%   samplestep is sampling step
%   projection is type of projection
% OUTPUT:
%   Acoarse is the resampling map where all "faces" are put row-wise
%   m is number of rectangular faces
%   n is the linear resolution of the faces
%   Acoarse is dimension [n] x [n*m] x [3]

switch projection
    case 'cube',
        % The topology of six faces is like this
        %        N
        %        ^
        %    W <-o-> E
        %        v
        %        S
        %      +---+
        %      | 6 |
        %  +---+---+---+---+
        %  | 4 | 1 | 2 | 3 |
        %  +---+---+---+---+
        %      | 5 |
        %      +---+
        %
        mx = 4;
        my = 3;
        
        subreg = [5 8 11 2 6 4];
    case 'mercator',
        % The topology of the face is like this (!)
        %        N
        %        ^
        %    W <-o-> E
        %        v
        %        S
        %      +---+
        %      | 1 |
        %      +---+
        %
        mx = 1;
        my = 1;
        subreg = [1]; %#ok
end

if isnan(nxy)
    nxy = [size(A,1) size(A,2)];
elseif isscalar(nxy)
    nxy = [nxy nxy];
end
% lower resolution
n = ceil(nxy./samplestep);
m = length(subreg);
j1 = (0:nxy(1)-1).';
j2 = (0:nxy(2)-1).';
if any(n.*samplestep > nxy) % this 'if' is not really needed
    % not exactly divided, pad with last index
    displayfun('Warning: resampling padding -> loss accuracy\n');
    j1(end+1:n(1)*samplestep) = nxy(1)-1;
    j2(end+1:n(2)*samplestep) = nxy(2)-1;
end
iy = bsxfun(@plus,j1,y1+((1:my)-1)*nxy(1)); % three faces on y direction
ix = bsxfun(@plus,j2,x1+((1:mx)-1)*nxy(2)); % four in x

% get submatrix
A = A(iy,ix,:);

% Get six-face only
nxy = n.*samplestep; % bug corrected
A = reshape(A, [nxy(1) my nxy(2) mx 3]);
A = reshape(permute(A, [1 3 2 4 5]), [nxy(1:2) mx*my 3]);
A = A(:,:,subreg,:); % six faces (nxy x nxy x 6 x 3)

szA = [n(1) n(2)*m];
% partitionning
A = reshape(A, [samplestep szA(1) samplestep szA(2) 3]);
% Put the subpixels in the last dimension
A = reshape(permute(A, [2 4 5 1 3]), szA(1), szA(2), 3, []);
% Cast and normalize to [0,1]
if strcmp(class(A),'uint8') % 0-255 RGB 
    A = double(A) / double(intmax(class(A)));
end
% Average the color on the new patch
szA = size(A);
Acoarse = zeros(szA(1:3),'double');
% for-loop is less memory demanding than one-shot mean
for k=1:size(A,3)
    Acoarse(:,:,k) = mean(A(:,:,k,:), 4);
end

end % resampling

%%
function IDX = GridIdx(Faces, FaceID, ID, n, projection)
% function IDX = GridIdx(Faces, FaceID, ID)
% Return grid index of one face (in a 2D array format)

switch projection
    case 'cube',
        % Filter
        ID = feval(class(FaceID),ID);
        F = Faces(FaceID==ID,:);
    case 'mercator',
        F = Faces; % Just one face!
end

F = reshape(F,[n(1) n(2) 4]);
IDX = [F(:,:,1)   F(:,end,2);
       F(end,:,4) F(end,end,3)];
   
end % GridIdx
