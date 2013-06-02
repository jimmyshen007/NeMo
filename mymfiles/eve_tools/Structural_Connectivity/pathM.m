function [path,a] = pathM(seed,seed_r_zone,step,lengthpath,factor_a,dir,interp)
% Usage : [path,a] = pathM(rAtlas,EV,Mask,seed,step,lengthpath,factor_a,dir,inf_Data)
%
% Description: This function compute nervous fibers trajectories front a set
% of seed points to a set of target gray matter zones using Mori's method.
%
% Input Parameters:
%        Atlas     : three-dimensional matrix that presents an specific value
%                     in each gray matter zone position.
%        Diff_Data         : diffusion data: four dimensional matrix whose first three
%                     dimensions indicate position and the fourth possess
%                     the principal eigenvector in each position.
%        Mask       : three-dimensional matrix that expresses the presence
%                     of white or gray matter in each point.
%        seed       : coordinates of the seed point.
%        step       : step size.
%        lengthpath : maximum number of steps that can have each path.
%        factor_a   : column vector, it depends to dimensions of the voxel
%                     corresponding to the diffusion images.
%        dir        :
%
% Output Parameters:
%        path       : path calculated.
%        a          : gray matter zones of which leave and arrive the path.
%
%--------------------------------------------------------------------------
% Useful References:
%  1.-Mori S., Crain J. B., Van Zijl P. C. M. y Chackov V. P. Three Dimensional
% Tracking of Axonal Projections in the Brain by Magnetic Resonance Imaging.
% Ann. Neurol, 45, 263-269, 1999.
%--------------------------------------------------------------------------
% Authors: Yasser Iturria Medina & Pedro Vald�s-Hern�ndez.
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

global Diff_Data inf_Data Mask Atlas

% path(:,1) = seed' + 0.5;
path(:,1) = seed';

%seednumber = Atlas(floor(seed(1)),floor(seed(2)),floor(seed(3)));
seednumber = seed_r_zone;
for stepn = 1:lengthpath - 1
    if interp==1 %%%%%%%%%%%%%%%%%Interp Tensor
        point = round(path(:,stepn));
        x = (point(1)-2:point(1)+2)'; y = (point(2)-2:point(2)+2)'; z = (point(3)-2:point(3)+2)';
        %         ind = find(x <= 0); x(ind) = []; ind = find(x > size(Mask,1)); x(ind) = [];
        %         ind = find(y <= 0); y(ind) = []; ind = find(y > size(Mask,2)); y(ind) = [];
        %         ind = find(z <= 0); z(ind) = []; ind = find(z > size(Mask,3)); z(ind) = [];
        [X Y Z] = meshgrid(x,y,z); X = X(:); Y = Y(:); Z = Z(:);
        V = [X Y Z];
        cant = size(x,1)*size(y,1)*size(z,1);
        for i = 1:cant
            if (V(i,1)<= 0) | (V(i,1)>size(Mask,1)) | (V(i,2)<= 0) | (V(i,2)> size(Mask,2))...
                    | (V(i,3)<= 0) | (V(i,3)> size(Mask,3))
                Dneib(:,:,i) = zeros(3,3);
            else
                Dneib(:,:,i) = inf_Data(:,:,V(i,1),V(i,2),V(i,3));
            end
        end
        Dp = zeros(3,3);
        pointR = path(:,stepn);
        Data = reshape(Dneib,[3 3 size(x,1),size(y,1),size(z,1)]);
        Dneib = [];
        Dp(1,1) = interp3(x,y,z,squeeze(Data(1,1,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(2,2) = interp3(x,y,z,squeeze(Data(2,2,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(3,3) = interp3(x,y,z,squeeze(Data(3,3,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(1,2) = interp3(x,y,z,squeeze(Data(1,2,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(1,3) = interp3(x,y,z,squeeze(Data(1,3,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(2,3) = interp3(x,y,z,squeeze(Data(2,3,:,:,:)),pointR(1),pointR(2),pointR(3),'spline');
        Dp(2,1) = Dp(1,2);
        Dp(3,1) = Dp(1,3);
        Dp(3,2) = Dp(2,3);
        [V L] = eig(Dp);
        e1 = V(:,3);
    else
        e1 = squeeze(Diff_Data(floor(path(1,stepn)),floor(path(2,stepn)),floor(path(3,stepn)),:));
        if stepn == 1, e1 = e1*dir; end
    end
    if isnan(sum(e1))
        path = []; a = [];
        return
    end
    if stepn > 1, e1 = e1*sign(e_ant'*e1); end
    e1 = e1/(norm(e1)+eps);
    e_ant = e1;
    next_point = path(:,stepn) + (step*factor_a.*e1);
    if (next_point(1) > size(Mask,1)) | (next_point(1) < 1) | ...
            (next_point(2) > size(Mask,2)) | (next_point(2) < 1) | ...
            (next_point(3) > size(Mask,3)) | (next_point(3) < 1)
        path = []; a = [];
        return
    elseif Mask(floor(next_point(1)),floor(next_point(2)),floor(next_point(3))) == 0
        path = []; a = [];
        return
    end
    path(:,stepn + 1) = next_point;
    endnumber = Atlas(floor(next_point(1)),floor(next_point(2)),floor(next_point(3)));
    %     if (endnumber ~= 0) & (endnumber ~= seednumber), break; end
    if (endnumber ~= seednumber), break; end %prueba
end
% if (endnumber ~= 0) & (endnumber ~= seednumber)
if (endnumber ~= seednumber) %prueba
    a = [seednumber endnumber];
else
    path = [];
    a = [];
end