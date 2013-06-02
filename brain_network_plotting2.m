function brain_network_plotting2(at, C, V, lobes, figstr,CColorMat,savefig_flg,savemovie_flg,maxnodesz,maxedgesz,sbi)
%inputs:
% at        3-D brain volume with labeled regions 1 to N (_Altas.img)
% C         N by N connectivity matrix with weights of connectivity
%           strength
% V         N by 1 vector of the node sizes. If you want uniform, just do
%           all ones.
% lobes     N by 1 vector with the lobe information.
% figstr    a string with the figure name to be saved
% CColorMat N by N matrix with 
% savefig_flg 0 or 1 depending on if you want to save the figures as .tif
% savemovie_flg 0 or 1 depending on if you want to save the movie.

dbstop if error
if nargin < 11 || isempty(sbi)
    sbi = [1 1 1 1];
end
if nargin < 10 || isempty(maxedgesz)
    maxedgesz = 0.5;
end
if nargin < 9 || isempty(maxnodesz)
    maxnodesz = 2.0;
end
if nargin < 8 || isempty(savemovie_flg)
    savemovie_flg = 1;
end
if nargin < 7 || isempty(savefig_flg)
    savefig_flg = 1;
end
if nargin < 6 || isempty(CColorMat)
    CColorMat = 55*ones(size(C));
end
if nargin < 5 || isempty(figstr)
    figstr = 'figure';
end
if nargin < 4 || isempty(lobes)
    lobes = ones(length(C),1);
end
if nargin < 3
    V = ones(length(C),1);
end

m = size(V,1);
if size(V, 2) > 1
    rgbflag = 1;
     nodelabels = max(abs(V), [], 2);    
    cmap = abs(V(:,1:3))./repmat(max(abs(V(:,1:3))), [m,1]);
else 
    rgbflag = 0;
    nodelabels = V;
end

%Can adjust the maximum size of the nodes and edges here.
sznode = maxnodesz/max(abs(nodelabels))*abs(nodelabels);
line_fact = maxedgesz/max(C(:)); 
figure(sbi(1));
subplot(sbi(2),sbi(3),sbi(4));
%figure;
hold on
axis equal

at(isnan(at)) = 0;
ss = size(at);
WB = ~(at(:)==0);
WB = double(reshape(WB,ss));
myvol_render(WB,[1 1 1],colormap(bone(5)),3)
alpha(.08)
lighting phong;

MAP=colormap(hsv(100));
GRAY = colormap(gray(10));

points=zeros(size(C,2),3);
kkc = 0;

for i = 1:size(C,2)  
    p  = [];
    [p(:,1),p(:,2),p(:,3)] = ind2sub(size(at),find(at(:)==i));    
    p = unique(p,'rows');
    points(i,:) = [mean(p(:,2)),mean(p(:,1)),mean(p(:,3))];
     % marker color by lobe name
     if length(unique(lobes))==6
         colors = {'blue';  'magenta'; 'green'; 'red'; 'cyan'; 'yellow'};
     elseif length(unique(lobes))==8
          colors = {'blue';  'magenta'; 'green'; 'red'; 'cyan'; 'yellow';'black';'white'};
     elseif length(unique(lobes))==2
          colors = {'blue';  'red'};
     elseif ~(length(unique(lobes))==6) && length(unique(lobes))>1
         nMAP = colormap(jet(240));
         nl = length(unique(lobes))-1;
         dvM = floor(size(nMAP,1)./nl);
         colors = nMAP(1:dvM:size(nMAP,1),:);
     end
    if nargin >=4 && (length(unique(lobes))==6 || length(unique(lobes))==8 || length(unique(lobes))==2)
        markercolor = colors{lobes(i)};
        if strcmp(markercolor,'white')
            markercolor = [1 0.5 0];
        end
    elseif nargin >=4 && ~(length(unique(lobes))==6) && length(unique(lobes))>1
        markercolor = colors(lobes(i),:);
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
            if C(j,i)>0
                kkc = kkc + 1;
                 P1 = [points(i,1),points(i,2), points(i,3)];
                 P2 = [points(j,1),points(j,2), points(j,3)];
                [XX,YY,ZZ] = cylinder2P(line_fact*C(i,j),100,P1,P2);
                if ~(CColorMat(j,i)==0)
                    surf(XX,YY,ZZ,'FaceColor',MAP(CColorMat(j,i),:),'EdgeColor','none')
                else
                    plot3([P1(1) P2(1)],[P1(2) P2(2)],[P1(3) P2(3)],'Color',GRAY(3,:),'LineWidth',line_fact*C(i,j))
                    %surf(XX,YY,ZZ,'FaceColor','k','EdgeColor','none')
                end
            end
        end
    end
end
fig = gcf;
set(fig, 'Name', figstr);  
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
set(gca,'Ztick',[]);
%xlim([16 112])
%ylim([27 106])

view([-90 90]);
camlight left
    
if savefig_flg
    % create views, videos and save them
    view([0 0 1]);
    camlight left
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_axial.tif']);
    clmo(handlem('light'))
    view([1 0 0]);
    camlight left
    pause(1);
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_coronal.tif']);
    clmo(handlem('light'))
    view([0 1 0]);
    camlight left
    pause(1);
    I = getframe(gcf);
    imwrite(I.cdata, [figstr '_sagittal.tif']);
    pause(1);
    clmo(handlem('light'))
end

if savemovie_flg
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
    movie2avi(M,[figstr '.avi'],'fps', 5); 
else
    camlight left
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
    %surf(X, Y, Z);