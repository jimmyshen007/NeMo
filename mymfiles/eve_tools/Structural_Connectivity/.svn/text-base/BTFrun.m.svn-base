load(['/home/evl2001/Documents/MATLAB/BTF/larabtf']);

diff_omit='';
diff_ext='s*.img';
% preatlasedoption=0;
atlassize=116;
bvecext='s*.bvec';
skipped=[];
preatlasedopt=0;
affinereg=1;

% for i=3:7
for i=8:14
    [~,f,~]=fileparts(deblank(larabtf(i,:)));
    PatDir=['/home/eve/Documents/MATLAB/BTFdata' filesep f(1:7)];
    PatDir=deb(larabtf(i,:));
    diff_dir=[PatDir filesep 'DTI'];
    T1FileName=deb(Dir2Arr([PatDir filesep 'T1scan'],'cos*.img'))
    if ~isempty(T1FileName)
        CalcPaths_Bayesian_SVN(PatDir,T1FileName,atlassize,diff_dir,diff_ext,diff_omit,bvecext,preatlasedopt,affinereg)
%         CalcPaths_Bayesian_SVN(PatDir,T1FileName,116,diff_dir,diff_ext,diff_omit,'*.bvec',1);
%         img_auto_translate(T1FileName,'backup_',116);
    else
        skipped=[skipped; f(1:7)];
    end
end

for i=1:21
    [~,f,~]=fileparts(deblank(evebtf(i,:)));
    PatDir=['/home/eve/Documents/MATLAB/BTFdata' filesep f(1:7)];
    tT1=deblank(Dir2Arr([PatDir filesep 'T1scan' filesep 'Normalized'],'w*.img'));
    openMRIcron(tT1,'t1-canonical');
end