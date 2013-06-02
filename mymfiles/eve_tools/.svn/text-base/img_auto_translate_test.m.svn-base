function img_auto_translate(coT1,RunAtTF,emlopt)

%The array coT1 points to the images that require shifting dcm2nii converted images to match T1
%template

%coT1=Dir2Arr([start_dir filesep 'test' filesep 'ANALYZE_translate_test'],'co*.img');

if nargin<2
    RunAtTF=0;
    emlopt=0;
end

%Set RunAtlasing command

RunAtCmd='RunAtlasing_307(face_off)';

if ischar(coT1), coT1=cellstr(coT1); end

tsp=[-0.5; -0.5872; -0.4066];

%T=size(coT1,1);
T=1;

try
for t=1:T
    [d,f,~]=fileparts(coT1{t});    
    V=spm_vol(coT1{t});
    for i=1:3
        V.mat(i,4)=tsp(i)*V.dim(i)*V.mat(i,i);
    end
    if ~exist([d filesep 'backup_' f '.hdr'],'file')
        copyfile([d filesep f '.hdr'],[d filesep 'backup_' f '.hdr']);
	copyfile([d filesep f '.img'],[d filesep 'backup_' f '.img']);
    end
    spm_create_vol(V);
    disp(f);
    face_off=Dir2Arr(d, ['*' f '.img']);
    if RunAtTF
        eval(RunAtCmd);
    end
end
   
catch
    if emlopt
        email_status(0,'img_auto_translate',['Job failed at t=' t]);
    end
end

if emlopt
    if t==T
        email_status(1,'img_auto_translate.m',['Job completed successfully.']);
    end
end
end
