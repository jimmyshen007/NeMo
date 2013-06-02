function Def = spm_calc_def_SVN(prm,flags)
% Write Out Warped Images.
% FORMAT VO = spm_write_sn(V,matname,flags,msk)
% V         - Images to transform (filenames or volume structure).
% matname   - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           interp   - interpolation method (0-7)
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
%           preserve - either 0 or 1.  A value of 1 will "modulate"
%                      the spatially normalised images so that total
%                      units are preserved, rather than just
%                      concentrations.
%           prefix   - Prefix for normalised images. Defaults to 'w'.
% msk       - An optional cell array for masking the spatially
%             normalised images (see below).
%
% Warped images are written prefixed by "w".
%
% Non-finite vox or bounding box suggests that values should be derived
% from the template image.
%
% Don't use interpolation methods greater than one for data containing
% NaNs.
% _______________________________________________________________________
%
% FORMAT msk = spm_write_sn(V,matname,flags,'mask')
% V         - Images to transform (filenames or volume structure).
% matname   - Transformation information (filename or structure).
% flags     - flags structure, with fields...
%           wrap     - wrap edges (e.g., [1 1 0] for 2D MRI sequences)
%           vox      - voxel sizes (3 element vector - in mm)
%                      Non-finite values mean use template vox.
%           bb       - bounding box (2x3 matrix - in mm)
%                      Non-finite values mean use template bb.
% msk       - a cell array for masking a series of spatially normalised
%             images.
%
%
% _______________________________________________________________________
%
% FORMAT VO = spm_write_sn(V,prm,'modulate')
% V         - Spatially normalised images to modulate (filenames or
%             volume structure).
% prm       - Transformation information (filename or structure).
%
%  After nonlinear spatial normalization, the relative volumes of some
%  brain structures will have decreased, whereas others will increase.
%  The resampling of the images preserves the concentration of pixel
%  units in the images, so the total counts from structures that have
%  reduced volumes after spatial normalization will be reduced by an
%  amount proportional to the volume reduction.
%
%  This routine rescales images after spatial normalization, so that
%  the total counts from any structure are preserved.  It was written
%  as an optional step in performing voxel based morphometry.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_sn.m 2534 2008-12-08 10:16:46Z christophe $

if ischar(prm), prm = load(prm);  end;

[x,y,z,mat] = get_xyzmat(prm,flags.bb,flags.vox);

[X,Y] = ndgrid(x,y);
Tr = prm.Tr;
BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1);
BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1);
BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1);

Def = [];
for j=1:length(z),   % Cycle over planes
    % Nonlinear deformations
    %----------------------------------------------------------------------------
    tx = get_2Dtrans(Tr(:,:,:,1),BZ,j);
    ty = get_2Dtrans(Tr(:,:,:,2),BZ,j);
    tz = get_2Dtrans(Tr(:,:,:,3),BZ,j);
    X1 = X    + BX*tx*BY';
    Y1 = Y    + BX*ty*BY';
    Z1 = z(j) + BX*tz*BY';
    
    [X2,Y2,Z2]  = mmult(X1,Y1,Z1,prm.Affine);
    Def.X(:,:,j) = X2;
    Def.Y(:,:,j) = Y2;
    Def.Z(:,:,j) = Z2;
end;
return;
%_______________________________________________________________________
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));
return;
%_______________________________________________________________________
function [X2,Y2,Z2] = mmult(X1,Y1,Z1,Mult)
if length(Z1) == 1,
    X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + (Mult(1,3)*Z1 + Mult(1,4));
    Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + (Mult(2,3)*Z1 + Mult(2,4));
    Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + (Mult(3,3)*Z1 + Mult(3,4));
else
    X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
    Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
    Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
end;
return;
%_______________________________________________________________________
function [x,y,z,mat] = get_xyzmat(prm,bb,vox,VG)
% The old voxel size and origin notation is used here.
% This requires that the position and orientation
% of the template is transverse.  It would not be
% straitforward to account for templates that are
% in different orientations because the basis functions
% would no longer be seperable.  The seperable basis
% functions mean that computing the deformation field
% from the parameters is much faster.

% bb  = sort(bb);
% vox = abs(vox);

if nargin<4,
    VG = prm.VG(1);

    if all(~isfinite(bb(:))) && all(~isfinite(vox(:))),
        x   = 1:VG.dim(1);
        y   = 1:VG.dim(2);
        z   = 1:VG.dim(3);
        mat = VG.mat;
        return;
    end
end

[bb0,vox0] = bbvox_from_V(VG);
if ~all(isfinite(vox(:))), vox = vox0; end;
if ~all(isfinite(bb(:))),  bb  = bb0;  end;

msk       = find(vox<0);
bb        = sort(bb);
bb(:,msk) = flipud(bb(:,msk));

% Adjust bounding box slightly - so it rounds to closest voxel.
% Comment out if not needed.
%bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
%bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
%bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = prm.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end;
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Convert range into range of voxels within template image
x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

og  = -vxg.*ogn;

% Again, chose whether to round to closest voxel.
%of  = -vox.*(round(-bb(1,:)./vox)+1);
of = bb(1,:)-vox;

M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = prm.VG(1).mat*inv(M1)*M2;

LEFTHANDED = true;
if (LEFTHANDED && det(mat(1:3,1:3))>0) || (~LEFTHANDED && det(mat(1:3,1:3))<0),
    Flp = [-1 0 0 (length(x)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
    mat = mat*Flp;
    x   = flipud(x(:))';
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [bb,vx] = bbvox_from_V(V)
vx = sqrt(sum(V.mat(1:3,1:3).^2));
if det(V.mat(1:3,1:3))<0, vx(1) = -vx(1); end;

o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
return;
%_______________________________________________________________________
