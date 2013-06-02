function openMRIcron(MRIimg, MNItemptype)

startup_varsonly;

if nargin<2
    MNItemptype=0;
end

if MNItemptype==0
    for i=1:size(MRIimg,1)
        system([mricronpath filesep 'mricron ' deblank(MRIimg(i,:)) ' &']);
    end
else
    for i=1:size(MRIimg,1)
        system([mricronpath filesep 'mricron ' deblank(MRIimg(i,:)) ' -o ' Choosing_TempFile(MNItemptype) ' &']);
    end
end

end