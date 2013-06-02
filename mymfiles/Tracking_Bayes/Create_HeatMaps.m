function Create_HeatMaps(qbelong_dir,ImageSpace)
% This function takes in the directory qbelong and outputs heat maps that
% are image volumes containing the number of fibers going through each voxel 
% that connect the two regions of interest. The mat file that is produced
% contains the indices of the non-zero voxels in the Heat map. For each
% nonzero voxel, there are two lists corresponding to the beginning and end 
% voxel for each of the N fibers in that voxel.

info = spm_vol(ImageSpace);
[~,~,ext] = fileparts(ImageSpace);
files = dir(qbelong_dir);
numROI = roots([1 -1 -2*(length(files)-2)]);
numROI = floor(numROI(find(numROI>0),1));

flag = strcmp(qbelong_dir(end-4:end),'NLMNI');
qbelongnew_dir = [qbelong_dir filesep 'HeatMaps'];
if ~isdir(qbelongnew_dir)
    mkdir(qbelongnew_dir)
end
tic

image_dim = info.dim;

for ii = 1:(numROI-1)
    disp(['starting region = ' num2str(ii)])
    for jj = ii+1:numROI
        load([qbelong_dir filesep 'qbelong' num2str(ii) '_' num2str(jj)]);
        heat_map = zeros(image_dim,'int32');
        if length(qbelong)==1
        else
            for j = 1:size(qbelong,2);
                p = qbelong(j).voxels;
                if flag;
                    p = p';
                end
                %pind = [];
                %Collect the voxels in the fibers.
                pind = sub2ind(info.dim,floor(p(1,2:end-1))',floor(p(2,2:end-1))',floor(p(3,2:end-1))');
                heat_map(pind) = heat_map(pind) + 1;
               % newcol_beg = NaN(prod(info.dim),1);
               % newcol_beg(pind) = sub2ind(info.dim,floor(p(1,1))',floor(p(2,1))',floor(p(3,1))');
               % pbeg = [pbeg; newcol_beg];
            end
            %Count the voxels and the number of fibers through them.
            hm_ind = find(heat_map); 
            num_fibers = nonzeros(heat_map);
            save([qbelongnew_dir filesep 'HeatMap' num2str(ii) '_' num2str(jj) '.mat'],'hm_ind','num_fibers','image_dim');
            %spm_create_vol(vinfo);
            %spm_write_vol(vinfo,heat_map);            
        end
    end
end