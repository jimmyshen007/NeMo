function [ACS, ACD, ACP] = Areas_saveConn(ROI1,ROI2,q,ROI_num1,ROI_num2,OutConPath)
% Usage : [ACS, ACD, ACP] = Area(ROI1,ROI2,q)
%
% Description: This function is to quantify anatomical connection strength
%              (ACS), anatomical connection density (ACD) and anatomical 
%              connection probability (ACP) between n gray matter zones using
%              information of computed routes that simulate nervous fibers 
%              trajectories.
%
% Input Parameters:
%        ROI1     : vector [n1,3], coordinates of the external voxels of the
%                   first gray matter zone.
%        ROI2     : vector [n2,3], coordinates of the external voxels of the
%                   second gray matter zone.
%        q        : variable struct that in the field 'voxels' possesses the 
%                   coordinates of each path that unite the two zones and in
%                   the field 'PT' possesses the probability or validity index
%                   of each of these paths. 
%  
% Output Parameters:
%        ACS      : Anatomical Connection Strength between ROI1 and ROI2.
%        ACD      : Anatomical Connection Density between ROI1 and ROI2.
%        ACP      : Anatomical Connection Probability between ROI1 and ROI2.
%--------------------------------------------------------------------------
% Useful References:
% 1.- Y. Iturria-Medina, P. Valdes-Hernandez, E. Canales-Rodriguez. 2005. 
% Measures of anatomical connectivity obtained from neuro diffusion images.
% Presented at the 11th Annual Meeting of the Organization for Human Brain
% Mapping, June 12-16, 2005, Toronto, Ontario, Canada. Available on CD-Rom
% in NeuroImage, Vol. 26, No.1.
%--------------------------------------------------------------------------
% Author: Yasser Iturria Medina & Pedro Vald�s-Hern�ndez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

if (~isempty(q(1).voxels)) && (~isempty(ROI1)) && (~isempty(ROI2))
    conn_z1_z2 = 0;
    conn_z2_z1 = 0;
    Ptos_zone1 = [];
    Ptos_zone2 = [];
    PT = [];
    Nf = size(q,2);
    for fibre = 1:Nf
        %disp(fibre)
        pt1 = floor(q(fibre).voxels(:,1))';
        pt2 = (q(fibre).voxels(:,end))';
        sdist = sum(abs(repmat(pt1,[size(ROI1,1) 1])-ROI1),2);
        i1 = min(sdist);
        sdist = sum(abs(repmat(pt1,[size(ROI2,1) 1])-ROI2),2);
        i2 = min(sdist);
        if i1 <= i2
            Ptos_zone1 = [Ptos_zone1; pt1];
            sdist = sum(abs(repmat(pt2,[size(ROI2,1) 1])-ROI2),2);
            [i,j] = min(sdist);
            Ptos_zone2 = [Ptos_zone2; ROI2(j,:)];
            PT(fibre) = abs(q(fibre).PT);
        elseif i1 > i2
            Ptos_zone2 = [Ptos_zone2; pt1];
            sdist = sum(abs(repmat(pt2,[size(ROI1,1) 1])-ROI1),2);
            [i,j] = min(sdist);
            Ptos_zone1 = [Ptos_zone1; ROI1(j,:)];
            PT(fibre) = abs(q(fibre).PT);
        end
    end
    if ~isempty(PT)
        Ptos_zone1_u = unique(Ptos_zone1,'rows');
        Ptos_zone2_u = unique(Ptos_zone2,'rows');
        for pt = 1:size(Ptos_zone1_u,1)
            ind = find((sum((repmat(Ptos_zone1_u(pt,:),[size(Ptos_zone1,1) 1])-Ptos_zone1),2))==0);
            conn_z1_z2 = conn_z1_z2 + max(PT(ind));
        end
        for pt = 1:size(Ptos_zone2_u,1)
            ind = find((sum((repmat(Ptos_zone2_u(pt,:),[size(Ptos_zone2,1) 1])-Ptos_zone2),2))==0);
            conn_z2_z1 = conn_z2_z1 + max(PT(ind));
        end
    end
    ACS = conn_z1_z2 + conn_z2_z1;
    ACD = conn_z1_z2/size(ROI1,1) + conn_z2_z1/size(ROI2,1);
    ACP = max([0 PT]);
    save([OutConPath filesep 'Ptos_' num2str(ROI_num1) '_' num2str(ROI_num2)],'PT','Ptos_zone1','Ptos_zone2')
else
    ACS = 0;
    ACD = 0;
    ACP = 0;
end