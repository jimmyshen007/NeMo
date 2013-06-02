function savePath=qbelong_to_trk(start_ind,qbdir,isofilename,atlassize)

startup_varsonly; 

if nargin<4
    atlassize=307;
end

isovol=spm_vol(isofilename);
modality=isovol.mat(1,1);


% [d,f,~]=fileparts(deblank(isofilename));
[d,f,~]=fileparts(deblank(isofilename));
[savedir,~,~]=fileparts(d);

savePath=[savedir filesep 'tracks_' start_ind '_' f '.trk'];

if ~exist(savePath,'file')
    disp(savePath);
    skiplen=length('qbelong')+length(start_ind)+1;
    idx=str2num(start_ind);
    %qblist=Dir2Arr(pbdir,{['qbelong' start_ind '_*.mat']; ['qbelong*_' start_ind '.mat']});
    qbstruct=[];
    qbnz=[];

    for i=1:(idx-1)
        load([qbdir filesep 'qbelong' int2str(i) '_' start_ind '.mat']);
        qbstruct(i).qb=qbelong;
        G=qbelong.voxels;
        if ~isempty(G)
            qbnz=[qbnz; i];
        end
    end

    for i=(idx+1):atlassize
        load([qbdir filesep 'qbelong' start_ind '_' int2str(i) '.mat']);
        qbstruct(i).qb=qbelong;
        G=qbelong.voxels;
        if ~isempty(G)
            qbnz=[qbnz; i];
        end
    end

    tracks=[];

    for j=1:length(qbnz)
        Q=qbstruct(qbnz(j)).qb.voxels;
        tracks(j).matrix=modality*Q';
%         tracks(j).matrix(:,1)=tracks(j).matrix(:,1)-20;
%         tracks(j).matrix(:,2)=tracks(j).matrix(:,2)-15;
%         tracks(j).matrix(:,3)=tracks(j).matrix(:,3)-15;
        tracks(j).nPoints=size(Q,2);
    end

%
%
%Load demo header to fill with our information
    exDir=[eve_tools filesep 'johncolby-along-tract-stats-db62641/example'];
    subDir  = fullfile(exDir, 'subject1');
    trkPath = fullfile(subDir, 'CST_L.trk');

    [header, ~] = trk_read(trkPath);

    T1vol=spm_vol(isofilename);

    %Update header.n_count with number of tracks loaded from qbelong
    header.n_count=length(tracks);
    header.dim=T1vol.dim;
    
    %set header.voxel_order
%     if T1vol.mat(1,1)/abs(T1vol.mat(1,1))==1
%         vox_ord(1)='R';
%     else
%         vox_ord(1)='L';
%     end
%     
%     if T1vol.mat(2,2)/abs(T1vol.mat(2,2))==1
%         vox_ord(2)='A';
%     else
%         vox_ord(2)='P';
%     end
%     
%     if T1vol.mat(3,3)/abs(T1vol.mat(3,3))==1
%         vox_ord(3)='S';
%     else
%         vox_ord(3)='I';
%     end
    
%     vox_ord(4)=' ';

   % vox
    
    for i=1:3
        eval('header.voxel_size(i)=abs(T1vol.mat(i,i));');
    end

    header.invert_y=1;
    header.voxel_order='LPS ';
    header.pad2=header.voxel_order;
    
    header.vox_to_ras=eye(4);
    
% Save the result back out to a .trk file for visualization in TrackVis.

    trk_write(header, tracks, savePath);

else
    disp(['Trk already exists: ' savePath]);
end

return;

end