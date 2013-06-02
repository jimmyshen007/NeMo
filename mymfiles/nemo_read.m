function [ PathROIs ] = nemo_read(fname, globalsize)
%NEMO_READ Summary of this function goes here
%   fname: input file name
%   globalsize: an array formated as [x y z] to represent the size of overall the
%   3D volume, used to index the voxel
    fid = fopen(fname, 'r');
    total = fread(fid, 1, 'int32');
    PathROI = struct('voxels', [], 'regions', [], 'location', []);
    PathROIs = struct(PathROI);
    for i=1:total
        vcols = fread(fid, 1, 'int32');
        vind = fread(fid, [1, vcols], 'int32');
        [rowsub colsub pagsub] = ind2sub(globalsize, vind);
        PathROI.voxels = [rowsub; colsub; pagsub];
        PathROI.regions = fread(fid, [1, 2], 'int16');
        PathROI.location = fread(fid, [1, 2], 'int16');
        PathROIs(1, i) = PathROI;
    end
    fclose(fid);
end

