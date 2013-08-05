function BrainographyRender(handles, axisHandle, scaleOverride)

if length(handles) == 1
    disp('Invalid rendering struct format');
    return;
end

if nargin > 2
    handles(end).renderRes = scaleOverride;
end

if nargin < 2
    axisHandle = figure;
end


renderProps = handles(1:end-1);
settingsStruct = handles(end);


axes(axisHandle);
set(gca,'DataAspectRatio',[1 1 1]);
hold on;

renderingVols = length(renderProps);

for i=1:renderingVols
    renderCortex(renderProps(i), axisHandle, settingsStruct);
   
    if renderProps(i).nodes || renderProps(i).pipes
        points = findCentroids(renderProps(i).brain_at, settingsStruct.renderRes);
    end
    
    if renderProps(i).nodes
        renderNodes(renderProps(i), points, settingsStruct);
    end
    
    if renderProps(i).pipes
        renderPipes(renderProps(i), points, settingsStruct);
    end
    
end

if settingsStruct.saveImages || settingsStruct.saveMovie
    if isempty(settingsStruct.figstr)
        figstr = 'brainography';
    else
        figstr = settingsStruct.figstr;
    end
    fig = gcf;
    set(fig, 'Name', figstr);
end

if settingsStruct.saveImages
    % create views, videos and save them
    view([0 0 1]);
    camlight left
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_axial.tif']);
    clmo(handlem('light'))
    view([1 0 0]);
    camlight left
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_sagittal.tif']);
    clmo(handlem('light'))
    view([0 1 0]);
    camlight left
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_coronal.tif']);
    clmo(handlem('light'))
end

if settingsStruct.saveMovie
%     zoom(1.5);
    [az,el] = view([1 0 0]);
    camlight left
    rot = [0:5:360];
    set(fig, 'Color', [1 1 1]);
    rotate3d;
    axis off;
    camproj perspective;
    axis vis3d;
    for i = 1:length(rot)
        view(rot(i), el);
        clmo(handlem('light'))
        camlight left
        M(i) = getframe(gcf);
    end 
    movie2avi(M,[figstr '.avi'],'fps', 5); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function renderCortex(renderProps, axisHandle, settingsStruct)

brain_at_bb = renderProps.brain_at;
brain_at_bb(isnan(brain_at_bb)) = 0;

% reduce resolution
if settingsStruct.renderRes ~= 1
    brain_at_bb = reduceRes(brain_at_bb, settingsStruct.renderRes);
end

opacity = renderProps.opacity;
% disp([num2str(opacity)]);
% present array structure: Column 1 is ROI-number, Column 2 is regionvalue,
% Column 3 is the color info as was calculated before passing to render;
% see self-contained assignRGB function for details

if ~isempty(renderProps.regionvalues)  %If regionvalue information provided (as rv or lut)
    
    numberROI = cell2mat(renderProps.regionvalues(:,1));
    regionValues = cell2mat(renderProps.regionvalues(:,2));
    rawMap = cell2mat(renderProps.regionvalues(:,3:5));
    
else  % Set default values if struct is empty
    
    numberROI = unique(brain_at_bb);
    numberROI = numberROI(find(numberROI));
    regionValues = numberROI;
    rawMap = getRGBTriple(bone(150),min(regionValues),max(regionValues),regionValues);
    
end

if renderProps.singleColorFlag
    binaryVol = zeros(size(brain_at_bb));
    for i=1:size(numberROI,1)
        if ~isnan(regionValues(i,:))
            binaryVol(brain_at_bb == numberROI(i)) = 1;
        end
    end
    [ccVol,~,numparts]=colorcode_regions(binaryVol);
    for k=1:numparts
        WB = zeros(size(brain_at_bb));
        WB(ccVol==k) = 1;
        renderCell(WB,rawMap(i,:),opacity);
    end
else
    for i = 1:size(numberROI,1)
        if ~isnan(regionValues(i,:))
            binaryVol = zeros(size(brain_at_bb));
            binaryVol(brain_at_bb == numberROI(i)) = 1;
            [ccVol,~,numparts]=colorcode_regions(binaryVol);
            for k=1:numparts
                WB = zeros(size(brain_at_bb));
                WB(ccVol==k) = 1;
                renderCell(WB,rawMap(i,:),opacity);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function renderCell(segment,ColorRGB,opacity)
% segment is a binary 3D volume to render, ColorRGB is the RGB triplet
% (0:1)

Ds = smooth3(segment*100);
hiso = patch(isosurface(Ds,5),'FaceColor',ColorRGB,'EdgeColor','none','FaceAlpha',opacity);
reducepatch(hiso,0.05);
isonormals(Ds,hiso);
% alpha(hiso,opacity);
lighting phong;

% set(hiso,'handlevisibility','off');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function renderNodes(renderProps, points, settingsStruct)

% C = renderProps.connectivityMatrix;
at = renderProps.brain_at;
% 
% if settingsStruct.renderRes ~= 1
%     reduceRes(at,settingsStruct.renderRes);
% end

if ~isempty(renderProps.nodeProps)
    numberROI = cell2mat(renderProps.nodeProps(:,1));
    if ~isempty(cell2mat(renderProps.nodeProps(:,2)))
        nodeVals = cell2mat(renderProps.nodeProps(:,2));
    else
        nodeVals = ones(size(numberROI));
    end
    if ~isempty(cell2mat(renderProps.nodeProps(:,3)))
        lobes = cell2mat(renderProps.nodeProps(:,3));
    else
        lobes = ones(size(numberROI));
    end
else
    numberROI = unique(find(unique(at)));
    nodeVals = ones(size(numberROI));
end


if ~isempty(renderProps.nodeScale)
    renderScale = renderProps.nodeScale/settingsStruct.renderRes;
else
    renderScale = 1.5/settingsStruct.renderRes;
end
sznode = renderScale/max(abs(nodeVals))*abs(nodeVals);

if isempty(renderProps.nodeSchema)
    colors = [0 0 1;  1 0 1 ; 0 1 0 ; 1 0 0 ; 0 1 1 ; 1 1 0];
else
    colors = renderProps.nodeSchema;
end

switch renderProps.nodeStyle
    case 1
        opacity = 1.0;
        singleRender = 1;
    case 2
        opacity = 0.25;
        singleRender = 1;
    case 3
        opacity(1) = 1.0;
        opacity(2) = 0.5;
        singleRender = 0;
end

for i = 1:length(lobes) 
    % marker color by lobe name
    if lobes(i) ~= 0 && sznode(i) ~= 0
        
        markercolor = colors(lobes(i),:);
        edgecolor = 'none';
        
        if singleRender
            [x,y,z]=ellipsoid(points(i,1), points(i,2), points(i,3), sznode(i), sznode(i), sznode(i), 200);
            surf(x,y,z,'FaceColor',markercolor,'EdgeColor',edgecolor,'FaceAlpha',opacity);
        else
            sznode0 = 0.7*sznode(i);
            [x0,y0,z0]=ellipsoid(points(i,1), points(i,2), points(i,3), sznode0, sznode0, sznode0, 200);
            [x1,y1,z1]=ellipsoid(points(i,1), points(i,2), points(i,3), sznode(i), sznode(i), sznode(i), 200);
            surf(x0,y0,z0,'FaceColor',markercolor,'EdgeColor',edgecolor,'FaceAlpha',opacity(1));
            surf(x1,y1,z1,'FaceColor',markercolor,'EdgeColor',edgecolor,'FaceAlpha',opacity(2));
        end
            
    end
end

function renderPipes(renderProps, points, settingsStruct)

C = renderProps.connectivityMatrix;


if renderProps.pipeUniform
    line_fact = (renderProps.pipeScale/settingsStruct.renderRes)/2.5;
else
    line_fact = (renderProps.pipeScale/settingsStruct.renderRes)/max(C(:));
end


switch renderProps.pipeScheme
    case 1 %single color
        CColorMat = renderProps.pipeColorHyperCube(:,:,:,1);
    case 2 %color map
        CColorMat = renderProps.pipeColorHyperCube(:,:,:,2);
    case 3 %2-color
        CColorMat = renderProps.pipeColorHyperCube(:,:,:,3);
end

if renderProps.pipeStyle==2
    opacity = 0.25;
else
    opacity = 1;
end

for i=1:size(C,2)
    if i>1
        for j=(i-1):-1:1
            if C(j,i) > 0
                P1 = [points(i,1), points(i,2), points(i,3)];
                P2 = [points(j,1), points(j,2), points(j,3)];
                if renderProps.pipeUniform
                    [XX,YY,ZZ] = cylinder2P(line_fact, 100, P1, P2);
                else 
                    [XX,YY,ZZ] = cylinder2P(line_fact*C(i,j), 100, P1, P2);
                end
                surf(XX,YY,ZZ, 'FaceColor', squeeze(CColorMat(j,i,:))', 'EdgeColor', 'none', 'FaceAlpha', opacity);
            end
        end
    end
end

function points = findCentroids(at, renderRes)

if renderRes ~= 1
    at = reduceRes(at,renderRes);
end

numberROI = unique(at(find(at ~= 0)));
points=zeros(size(numberROI,1),3);

for i = 1:size(numberROI,1)  
    tmp = zeros(size(at));
    tmp(at == numberROI(i)) = 1;
    p = regionprops(tmp, 'Centroid');
    
    points(i,:) = [p.Centroid(1) p.Centroid(2) p.Centroid(3)];
%     points(i,:) = [p.Centroid(2) p.Centroid(1) p.Centroid(3)];
end

function [X, Y, Z] = cylinder2P(R, N,r1,r2)

% The parametric surface will consist of a series of N-sided
% polygons with successive radii given by the array R.
% Z increases in equal sized steps from 0 to 1.

% Set up an array of angles for the polygon.
theta = linspace(0,2*pi,N);

m = length(R);                 % Number of radius values
% supplied.

if m == 1                      % Only one radius value supplied.
    R = [R; R];                % Add a duplicate radius to make
    m = 2;                     % a cylinder.
end


X = zeros(m, N);             % Preallocate memory.
Y = zeros(m, N);
Z = zeros(m, N);

v=(r2-r1)/sqrt((r2-r1)*(r2-r1)');    %Normalized vector;
%cylinder axis described by: r(t)=r1+v*t for 0<t<1
R2=rand(1,3);              %linear independent vector (of v)
x2=v-R2/(R2*v');    %orthogonal vector to v
x2=x2/sqrt(x2*x2');     %orthonormal vector to v
x3=cross(v,x2);     %vector orthonormal to v and x2
x3=x3/sqrt(x3*x3');

r1x=r1(1);r1y=r1(2);r1z=r1(3);
r2x=r2(1);r2y=r2(2);r2z=r2(3);
vx=v(1);vy=v(2);vz=v(3);
x2x=x2(1);x2y=x2(2);x2z=x2(3);
x3x=x3(1);x3y=x3(2);x3z=x3(3);

time=linspace(0,1,m);
for j = 1 : m
    t=time(j);
    X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x;
    Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y;
    Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
end

function redAt = reduceRes(at,renderRes) 
redAt = at(1:renderRes:end,1:renderRes:end,1:renderRes:end);
