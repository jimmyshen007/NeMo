function [] = PlotChaCoResults_ORIG(ChaCoResultsFileName,GBPlot,SurfPlot,BoxPlot,GraphPlot,figstr,plotlobecolor)

[pth,fn,~] = fileparts(ChaCoResultsFileName);
load(ChaCoResultsFileName)
if strfind(fn,'_T1')
    strsave = ['T1_' figstr];
else
    strsave =figstr;
end
% load the atlas of a subject. 
atsz = length(ChaCoResults(1).Regions);
AtlasName = ['/home/amy/work/BTFdata/Normals/e007573/T1scan/Atlased' num2str(atsz) '/sacos007a1001_Atlas.img'];
at = spm_read_vols(spm_vol(AtlasName));
at(isnan(at)) = 0;

if atsz==116 && plotlobecolor
    % list of lobe regions
    % 1 = frontal, 2 = parietal, 3 = occipital, 4 = temporal, 5 = subcortical,
    % 6 = cerebellum
    lobes = zeros(116,1);
    lobes(1:28) = 1;
    lobes([29:40, 71:78]) = 5;
    lobes(43:54) = 3;
    lobes([55:56, 79:90]) = 4;
    lobes(57:70) = 2;
    lobes(91:116) = 6;
    % new mods in limbic system
    lobes([29, 30]) = 1;  % insula
    lobes([31, 32]) = 1;  % ant cing
    lobes([33, 34]) = 1; % mid cing
    lobes([35, 36]) = 2; % post cing
    lobes([37, 38]) = 5; % hippo
    lobes([39, 40]) = 4; % parahippo 
    lobes([41, 42]) = 4; % parahippo
    %load /home/amy/work/WMH_Glodzik/CoregMasks/MeanCCL_13pats.mat
    C = zeros(116,116);
else
    lobes = ones(86,1);
    
    C = zeros(86,86);
end
if plotlobecolor==0
    lobes(find(ChaCoResults(1).Regions>0)) = 3;
    lobes(find(ChaCoResults(1).Regions<0)) = 4;
end

ZonesName = ReadInTxt(['/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Atlasing/atlas' num2str(atsz) '.cod']);
 
%plot the glass brain
if GBPlot.flag ==1
    brain_network_ploting(pth,at, C, abs(ChaCoResults(1).Regions'), [], lobes, 'GBmeanChaCo',[],GBPlot.movie,strsave)
end
%plot the gummi brain
if SurfPlot.flag == 1
    PlotIndex = floor((atsz-1)*(ChaCoResults(1).Regions'./max(ChaCoResults(1).Regions'))) + 1;
    GummiBrain(pth,PlotIndex,AtlasName,SurfPlot.MAP,SurfPlot.PlotHemi,ZonesName,lobes,'SPmeanChaCo',SurfPlot.movie,strsave)
end
%plot the boxplots of the ChaCo scores across the normal subjects
if BoxPlot.flag == 1
    mL = [];
    if length(ChaCoResults)==1
        mL = ChaCoResults(1).AllRegions;
    else
        for i = 2:length(ChaCoResults);
            mL = [mL ChaCoResults(i).Regions'];
        end
    end
    %map = lines(length(unique(lobes)));
    %colors = {'blue';  'magenta'; 'green'; 'red'; 'cyan'; 'yellow'};
    map = [0 0 1; 1 0 1; 0 1 0; 1 0 0; 0 1 1; 1 1 0];
    bp = figure; 
    screenSize = get(0,'ScreenSize');
    gbp = get(bp,'Position');
    set(bp, 'Position', [1 gbp(2) screenSize(3)-10 floor(screenSize(4)/2)])
     cmapbox = map(lobes,:);
    h = []; p = [];
    for i = 1:size(mL,1);
        [p(i),h(i)] = signtest(mL(i,:),[],'alpha',0.05/size(mL,1),'method','exact');
    end
    cmapbox(find(h==0),:) = repmat([0 0 0],[length(find(h==0)) 1]);
    boxplot(mL','plotstyle','compact','colors',cmapbox,'labels',ZonesName) 
    I = getframe(bp);
    imwrite(I.cdata, [pth filesep 'ChaCoPerRegion' strsave '.tif']);
end
if GraphPlot.flag == 1
    mu_OM = load(['/mnt/extra/BTFdata/FiberTracts' num2str(atsz) '_MNI/AllConnMatrices' num2str(atsz) '_FT_MNI/nGraphMets' num2str(atsz) '_MNI']);
    netmet = [];
    netmet_NC = [];
    if length(ChaCoResults)==1
        netmet = ChaCoResults.AllnConMatMets;
        netmet_NC = ChaCoResults(1).OrigMat.AllnConMatMets;
    end    
    bp = figure; 
    screenSize = get(0,'ScreenSize');
    gbp = get(bp,'Position');
    set(bp, 'Position', [1 gbp(2) floor(screenSize(3)/2) floor(screenSize(4)/2)])
    subplot(1,2,1)
    title('Characteristic Path Length','FontSize',16)
    hold on
    boxplot(netmet(:,1),'colors','r','labels',{'Patient'})
    hold on;
    %plot(2,netmetm(1),'r.','MarkerSize',24)
    %hold on
    plot(1,netmet_NC(1),'b.','MarkerSize',24)
    hold on
    legend({'Mean of NC Matrices'},'Location','SouthWest','FontSize',12)
    hold on;
   % ylim([1525 max(netmet(:,1))+25])
    hold off;    
    subplot(1,2,2)
    title('Efficiency','FontSize',16)
    hold on
    boxplot(netmet(:,2),'colors','r','labels',{'Patient'})
    hold on;
    %plot(2,netmetm(2),'r.','MarkerSize',24)
    %hold on
    plot(1,netmet_NC(2),'b.','MarkerSize',24)
    hold on
    legend({'Mean of NC Matrices'},'Location','SouthWest','FontSize',12)
    hold on;
    ylim([min(8.5e-4,min(netmet(:,2))) 9.4e-4])
    hold off;
    I = getframe(gcf);
    imwrite(I.cdata, [pth filesep 'GraphMets_ConMat' strsave '.tif']);
end
%--------------------------------------------------------------------------
function brain_network_ploting(pth,at, C, V, roi_names, lobes, figstr,CColorMat,movie_flg,strsave)

dbstop if error
MAP=colormap(hsv(100));
m = size(V,1);
if size(V, 2) > 1
    rgbflag = 1;
 
     nodelabels = max(abs(V), [], 2);    
  
    cmap = abs(V(:,1:3))./repmat(max(abs(V(:,1:3))), [m,1]);
else 
    rgbflag = 0;
    nodelabels = V;
end
if nargin >6
else
    figstr = 'figure';
end
sznode = 4/max(abs(nodelabels))*abs(nodelabels);
line_fact = 1.0/max(C(:)); 
bp = figure('visible','off'); 
screenSize = get(0,'ScreenSize');
gbp = get(bp,'Position');
set(bp, 'Position', [floor(1*screenSize(3)/2) screenSize(4) floor(3*screenSize(4)/4) floor(1*screenSize(4)/2)])
hold on
axis equal

at(isnan(at)) = 0;
ss = size(at);
WB = ~(at(:)==0);
WB = double(reshape(WB,ss));
myvol_render(WB,[1 1 1],colormap(bone(5)),3)
alpha(.08)
lighting phong;

points=zeros(size(C,2),3);
kkc = 0;
k32 = 0;
k100 = 0;

for i = 1:size(C,2)  
    p  = [];
    [p(:,1),p(:,2),p(:,3)] = ind2sub(size(at),find(at(:)==i));    
    p = unique(p,'rows');
    points(i,:) = [mean(p(:,2)),mean(p(:,1)),mean(p(:,3))];
     % marker color by lobe name
    colors = {'blue';  'magenta'; 'green'; 'red'; 'cyan'; 'yellow'};
    if nargin >4
        markercolor = colors{lobes(i)};
    else
        markercolor = 'red';
    end
    if nodelabels(i) > 0
        edgecolor = 'none';
    else 
        edgecolor = 'none';
    end
    [x,y,z]=ellipsoid(mean(p(:,2)),mean(p(:,1)),mean(p(:,3)),sznode(i),sznode(i),sznode(i),200);
    if rgbflag ==1
        markercolor = cmap(i,:);
        if sznode(i) > 0
            surf(x,y,z,'FaceColor',markercolor,'EdgeColor',edgecolor)
       end     
    else
        if sznode(i) > 0
            surf(x,y,z,'FaceColor',markercolor,'EdgeColor',edgecolor)
       end
    end
    
    if i>1
        for j=(i-1):-1:1
            if C(j,i)>0%.03*max(C(:))
                kkc = kkc + 1;
                 P1 = [points(i,1),points(i,2), points(i,3)];
                 P2 = [points(j,1),points(j,2), points(j,3)];
                [XX,YY,ZZ] = cylinder2P(line_fact*C(i,j),100,P1,P2);
                % [XX1,YY1,ZZ1]= tubeplot(XX,YY,ZZ,line_fact*C(j,i),[],2);
                if CColorMat(j,i)==32; 
                    k32 = k32 + 1;
                elseif CColorMat(j,i) == 100
                    k100 = k100 + 1;                    
                end
                surf(XX,YY,ZZ,'FaceColor',MAP(CColorMat(j,i),:),'EdgeColor','none')
                %plot3(points([i,j],1),points([i,j],2), points([i,j],3), 'k-',  'LineWidth', (line_fact*C(j,i)));
            end
        end
    end
end
fig = gcf;
set(fig, 'Name', figstr);  

% create views, videos and save them
%     scrsz = get(0,'ScreenSize');
%     set(fig, 'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
view([0 0 1]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_axial' strsave '.tif']);
clmo(handlem('light'))
view([0 1 0]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_sagittal' strsave '.tif']);
clmo(handlem('light'))
view([1 0 0]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_coronal' strsave '.tif']);

if movie_flg
    zoom(2);
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
    movie2avi(M,[pth filesep figstr '_' strsave '.avi'],'fps', 5); 
    zoom('reset')
end

return

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
return

function [] = myvol_render(seg,side_len,Map,Color)
% seg is a binary 3D volume to render.  side_len is the aspect ratio of
% each dimension. that we want to render.
[h,w,d] = size(seg);
Ds = smooth3(seg*100);
hiso = patch(isosurface(Ds,5),'FaceColor',Map(Color,:),'EdgeColor','none');
isonormals(Ds,hiso);
view(35,30);
axis([1 w 1 h 1 d]);
side_len = side_len([2 1 3]);
daspect(1./side_len);
return

function GummiBrain(pth,ColorIndex,AtlasName,MAP,PlotHemi,ZonesName,lobes,figstr,movie_flag,strsave)
%ColorIndex is the intensity of the color you want to plot for each
%region, you need to provide this. It will have one number per region in
%your atlas.

%AtlasName is the atlased image that you want to plot the brain on 
%(this is the atlased T1 image in the "Atlased" subfolder of the T1scan 
%subdirectory).
at = spm_read_vols(spm_vol(AtlasName));
at(isnan(at)) = 0;
ss = size(at);
roinums = setdiff(unique(at(:)),0);

%MAP: the colormap you want to use. just check the help for colormap if you
%want to change this.
if nargin<3 || isempty(MAP)
    MAP = colormap(hot(length(roinums)+1));
    MAP = MAP(1:end-1,:);
end
%If you want to plot the whole brain, use the default setting or input
%PlotHemi as a vector of all ones; If you want to plot just the left or right hemisphere,
%input PlotHemi as a binary vector indicating the regions you want to plot
% with 1's.
bp = figure; 
screenSize = get(0,'ScreenSize');
gbp = get(bp,'Position');
set(bp, 'Position', [floor(1*screenSize(3)/2) 1 floor(3*screenSize(4)/4) floor(1*screenSize(4)/2)])
hold on
if nargin<4
    PlotHemi = 'both';
end
if strcmp(PlotHemi,'subcortical')
    at(isnan(at)) = 0;
    ss = size(at);
    WB = ~(at(:)==0);
    WB = double(reshape(WB,ss));
    myvol_render(WB,[1 1 1],colormap(bone(5)),3)
    alpha(.08)
    lighting phong;
    hold on;
end
for i = roinums';    
    if strcmp(PlotHemi,'left') && ~isempty(strfind(ZonesName(i,:),'_L')) && isempty(strfind(ZonesName(i,:),'_Lobule_R'))
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    elseif strcmp(PlotHemi,'right') && ~isempty(strfind(ZonesName(i,:),'_R'))
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    elseif strcmp(PlotHemi,'subcortical') && lobes(i)==5
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    elseif strcmp(PlotHemi,'both')
        plot3d = at(:)==i;
        plot3d = double(reshape(plot3d,ss));
        myvol_render(plot3d,[1 1 1],MAP,ColorIndex(i))
        hold on;
    end    
end
lightangle(45,30);
lighting phong;
figstr = [figstr '_' PlotHemi];

fig = gcf;
set(fig, 'Name', figstr);

% create views, videos and save them
%     scrsz = get(0,'ScreenSize');
%     set(fig, 'OuterPosition',[1 1 scrsz(3) scrsz(4)]);
view([0 0 1]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_axial' strsave '.tif']);
clmo(handlem('light'))
view([0 1 0]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_sagittal' strsave '.tif']);
clmo(handlem('light'))
view([1 0 0]);
camlight left
pause(0.1)
I = getframe(gcf);
imwrite(I.cdata, [pth filesep figstr '_coronal' strsave '.tif']);
if strcmp(PlotHemi,'left') 
    clmo(handlem('light'))
    view([0 -1 0]);
    camlight left
    pause(0.1)
    I = getframe(gcf);
    imwrite(I.cdata, [pth filesep figstr '_medial' strsave '.tif']);
end
if strcmp(PlotHemi,'right')
    clmo(handlem('light'))
    view([0 -1 0]);
    camlight left
    pause(0.1)
    I = getframe(gcf);
    imwrite(I.cdata, [pth filesep figstr '_sagittal' strsave '.tif']);
    clmo(handlem('light'))
    view([0 1 0]);
    camlight left
    pause(0.1)
    I = getframe(gcf);
    imwrite(I.cdata, [pth filesep figstr '_medial' strsave '.tif']);
end
%clmo(handlem('light'))

if movie_flag
    zoom(2);
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
    movie2avi(M,[pth filesep figstr strsave '.avi'],'fps', 5); 
end
return
