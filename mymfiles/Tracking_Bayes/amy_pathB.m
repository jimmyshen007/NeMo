function [path,a,Ptm_T] = amy_pathB(seed,seed_zones,lengthpath,factor_a,ang)
% Usage : [path,a,Ptm_T] = pathB(rAtlas,Diff_Data,Dir0,Mask,seed,lengthpath,factor_a,ang)
%
% Description: This function compute nervous fibers trajectories front a set
% of seed points to a set of target gray matter zones using a Bayesian
% Formulation ("BayesianF"). Also, it compute the probability or validity
% index of each path.
%
% Input Parameters:
%        Atlas      : three-dimensional matrix that presents an specific value
%                     in each gray matter zone position.
%        seed       : matrix [n,3] that contains the coordinates of the n seed
%                     points.
%        Diff_Data    : diffusion data: ODF or four dimensional matrix whose first
%                     three dimensions indicate position and the fourth possess
%                     the value of fiber distribution for each one of the M
%                     directions.
%        inf_Data       : complementary information on the diffusion data:
%                     [M,3] matrix that contains the M unitary vectors
%                     of directions of the ODF.
%        Mask       : three-dimensional matrix that expresses the presence
%                     of white or gray matter in each point.
%        seed       : coordinates of the seed point.
%        lengthpath : maximum number of steps that can have each path.
%        factor_a   : column vector, it depends to dimensions of the voxel
%                     corresponding to the diffusion images.
%        step       : step size.
%        ang        : maximum angle of bend that can have each path
%                     between two serial steps.
%
% Output Parameters:
%        path       : path calculated.
%        a          : gray matter zones of which leave and arrive the path.
%        Ptm_T      : probability of the path.
%
%--------------------------------------------------------------------------
% Useful References:
%  1.-Y. Iturria-Medina, E. Canales-Rodríguez, L. Melié-García, P. Valdés-Hernández.
% Bayesian formulation for fiber tracking. Presented at the 11th Annual Meeting
% of the Organization for Human Brain Mapping, June 12-16, 2005, Toronto, Ontario,
% Canada. Available on CD-Rom in NeuroImage, Vol. 26, No.1.
%--------------------------------------------------------------------------
% Authors: Yasser Iturria Medina & Erick Canales Rodríguez.
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

global Diff_Data inf_Data Mask Atlas

Ptm_T = 1;
path(:,1) = seed' + 0.5;
%seednumber = Atlas(floor(seed(1)),floor(seed(2)),floor(seed(3)));
seednumber = seed_zones;
PT_o = 0;
do = 0;f_o = 0;
stepn = 0; long = 0;
if ndims(Diff_Data) == 5, tensor=1; else tensor=0; end
n = size(inf_Data,1);
while long < lengthpath
    stepn = stepn + 1;
    point = floor(path(:,stepn));
    if tensor == 1
        DTensor = squeeze(Diff_Data(:,:,point(1),point(2),point(3)));
        if sum(sum(DTensor)) ~= 0
            Dinv = inv(DTensor);
            odfv = abs(diag(inf_Data*Dinv*inf_Data')).^(-0.5);
            odfv = odfv/(sum(odfv) + eps);
            %             GFA = ((n*sum((odfv - mean(odfv)).^2))/((n-1)*sum(odfv.^2))).^0.5;
            %             f = GFA*(odfv-min(odfv))/(max(odfv)-min(odfv));
            f = (odfv-min(odfv))/(max(odfv)-min(odfv));
            f = f/sum(f + eps);
        else
            path = []; a = []; Ptm_T = [];
            return
        end
    else
        f = squeeze(Diff_Data(point(1),point(2),point(3),:));
    end
    f0 = f;
    if sum(f) == 0
        path = []; a = []; Ptm_T = [];
        return
    end
    ind = find(f < 10^(-3));
    f(ind) = [];
    if stepn > 1
        c = 0;
        MatrixSig = sign(inf_Data * do);
        Dir = repmat(MatrixSig,[1 3]).*inf_Data;
        Dir(ind,:) = [];
        Cos = Dir * (do/norm(do));
        ind = find(Cos < cos(ang));
        Dir(ind,:) = [];
        f(ind) = [];
        Do = repmat(do',[size(f,1) 1]);
    else
        c = 1;
        Dir = inf_Data;
        Dir(ind,:) = [];
        Do = 0;
    end
    if isempty(f)
        path = []; a = []; Ptm_T = [];
        return
    end
    N = size(Dir,1);
    D = ((c+1)/2)*(Dir.*repmat(f/max(f),[1 3]) + (1-c)*PT_o*Do);
    Mat_norm = sqrt(diag(D*D'))';
    factor = f/max(f);
    q = repmat(path(:,stepn)',[N 1]) + (D./(repmat(Mat_norm',[1 3])+eps)).*repmat(factor,[1 3]).*repmat(factor_a',[N 1]);
    Q = floor(q);
    ind = find(Q(:,1) > size(Mask,1)); ind = [ind ; find(Q(:,1) < 1)];
    ind = [ind ;find(Q(:,2) > size(Mask,2))]; ind = [ind ; find(Q(:,2) < 1)];
    ind = [ind ;find(Q(:,3) > size(Mask,3))]; ind = [ind ; find(Q(:,3) < 1)];
    ind = unique(ind);
    %     if size(ind,1) == size(q,1)
    %         break
    %     end
    q(ind,:) = []; Q(ind,:) = [];
    f(ind) = []; Dir(ind,:) = [];
    Ptm = [];  D(ind,:) = [];
    for i = 1:size(Q,1)
        Ptm(i) = Mask(Q(i,1),Q(i,2),Q(i,3));
    end
    if sum(Ptm) == 0
        path = []; a = []; Ptm_T = [];
        return
    end
    %%%%%%%%%%
    [m,Posc] = max(abs(inf_Data*D'),[],1);
    f_D = f0(Posc);
    %%%%%%%%%%
    PT = Ptm.*f_D';
    if sum(PT) == 0
        path = []; a = []; Ptm_T = [];
        return
    end
    PT = PT/sum(PT);
    PT_suma = [];
    for i = 1:size(PT,2)
        PT_suma(i) = sum(PT(1,1:i));
    end
    a = rand;
    vector_diff = PT_suma - a;
    ind = min(vector_diff(find(vector_diff > 0)));
    ind = find(vector_diff == ind);
    path(:,stepn + 1) = q(ind(1),:)';
    if (path(1,stepn + 1) > size(Mask,1)) | (path(1,stepn + 1) < 1) | ...
            (path(2,stepn + 1) > size(Mask,2)) | (path(2,stepn + 1) < 1) | ...
            (path(3,stepn + 1) > size(Mask,3)) | (path(3,stepn + 1) < 1)
        path = []; a = []; Ptm_T = [];
        return
    elseif Mask(round(path(1,stepn + 1)),round(path(2,stepn + 1)),round(path(3,stepn + 1))) == 0
        path = []; a = []; Ptm_T = [];
        return
    end
    do = path(:,stepn + 1) - path(:,stepn);
    PT_o = PT(ind(1));
    %     Ptm_T = ((Ptm_T^(stepn-1))*PT_o)^(1/stepn);
    Ptm_T = min([Ptm_T PT_o]);
    endnumber = Atlas(Q(ind(1),1),Q(ind(1),2),Q(ind(1),3));
    if (endnumber ~= 0) & (endnumber ~= seednumber)
        break
    end
    long = long + norm(do);
end
if (endnumber ~= 0) & (endnumber ~= seednumber)
    a = [seednumber endnumber];
else
    path = [];
    a = [];
    Ptm_T = [];
end