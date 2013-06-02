function ratbatch(t1,atlassize,MRIOPT)

if nargin < 3
    MRIOPT=0;
    if nargin < 2
        atlassize=116;
    end
end

for i=1:size(t1,1)
    disp([num2str(i) ' out of ' num2str(size(t1,1)) ' images']);
    RunAtlasing_generic(t1(i,:),atlassize,MRIOPT);
end
end