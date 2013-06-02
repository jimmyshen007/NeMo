function avgmatrix=qbelong_trklength_avg_SVN(qbdir,atlassize)

if nargin<2
    atlassize=307;
end

avgmatrix=zeros(atlassize,atlassize);
   
for i=1:atlassize
    %disp(['i=' num2str(i)]);
    for j=i+1:atlassize
        %disp(['j=' num2str(j)]);
        load([qbdir filesep 'qbelong' num2str(i) '_' num2str(j) '.mat']);
        L=length(qbelong);
        net_legs=0;
        if ~isempty(qbelong(1).voxels)
            for l=1:L
            %disp(['l=' num2str(l)]);
                V=size(qbelong(l).voxels,2);
                trackdist=0;
                %Scale coordinates by voxel dimension (1.8 in this case)
                testtrack=1.8*qbelong(l).voxels';
                for v=2:V
                %disp(['v=' num2str(v)]);
                    leg=sqrt((testtrack(v-1,1)-testtrack(v,1))^2+(testtrack(v-1,2)-testtrack(v,2))^2+(testtrack(v-1,3)-testtrack(v,3))^2);
                    %leg=(sum(testtrack(v-1,:)-testtrack(v,:)).^2)^(1/2);
                    trackdist=trackdist+leg;
                end
                net_legs=net_legs+trackdist;
            end
            avgmatrix(i,j)=net_legs/L;
            avgmatrix(j,i)=avgmatrix(i,j);
        end
    end
end
[savedir,~,~]=fileparts(qbdir);
save([savedir filesep 'avg_path_lengths' num2str(atlassize) '.mat'],'avgmatrix');
figure; imagesc(avgmatrix);
return;
end
