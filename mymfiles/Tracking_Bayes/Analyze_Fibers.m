function [Map,L] = Analyze_Fibers(qbelong_dir,DiffImage)
%This function takes in the directory qbelong and outputs the number of
%fibers between two regions (contained in the N x N matrix C), the number
% of voxels with fibers passing through them between two regions and a
%vector of fiber lengths (in mm) of all the fibers in qbelong.

info = spm_vol(DiffImage);
vox_sz = info.mat(1,1);
files = dir(qbelong_dir);
numROI = roots([1 -1 -2*(length(files)-2)]);
numROI = floor(numROI(find(numROI>0),1));
C = zeros(numROI);
fmax = 10000;
L = NaN(numROI,numROI,fmax);
Lmean = zeros(numROI);
Lstd = zeros(numROI);
N = zeros(numROI);
Con = [];
tic
for ii = 1:numROI
    disp(['starting region = ' num2str(ii)])
    for jj = ii+1:numROI
        load([qbelong_dir filesep 'qbelong' num2str(ii) '_' num2str(jj)]);
        if isempty(qbelong(1).voxels)
            C(ii,jj) = 0;
        else
            C(ii,jj) = size(qbelong,2);
            if size(qbelong,2)>fmax; L(:,:,fmax+1:size(qbelong,2)) = NaN; fmax = size(qbelong,2); end
            leng = [];
            pind = [];
            for j = 1:size(qbelong,2);
                p = qbelong(j).voxels;
                %Collect the voxels in the fibers.
                pind = [pind; sub2ind(info.dim,floor(p(1,:))',floor(p(2,:))',floor(p(3,:))')];
                diff = NaN(size(p,2)-1,1);
                for k = 2:size(p,2)
                    diff(k-1) = sqrt(sum((p(:,k-1)-p(:,k)).^2));
                end
                leng(j) = vox_sz*sum(diff);
            end
            %Count the voxels and the number of fibers through them.
            [Con(ii,jj).vox,Con(ii,jj).count] = count_unique(pind);
            N(ii,jj) = length(unique(pind));
            L(ii,jj,1:size(leng,2)) = leng;
            if ~(size(leng,2)==size(qbelong,2)); disp('FUCK!!'); end
            Lmean(ii,jj) = mean(leng);
            Lstd(ii,jj) = std(leng);
        end
    end
end
time_to_count = toc
Map = C + C';
N_vox = N + N';
L = single(L);
save([qbelong_dir(1:end-8) filesep 'Matrix_FC.mat'],'Map','L','Lmean','Lstd','Con','N_vox')