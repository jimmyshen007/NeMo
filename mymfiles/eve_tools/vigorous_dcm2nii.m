function vigorous_dcm2nii(ImgSrc,inifile,ext)

startup;
%T2src='/home/eve/Documents/T2_DCM';

dirlist=dir(ImgSrc);
%inifile='/home/eve/Documents/MATLAB/t2_dcm2nii.ini';
logfile=[eve_tools filesep 'd2n_log.txt'];

niigz_list=[];

for i=3:length(dirlist)
    indir=[ImgSrc filesep dirlist(i).name];
    disp(['system call to ''dcm2nii'': Converting ' indir ' ...']);
    system(['dcm2nii -b ' inifile ' ' indir ' > ' logfile]);
    
    logfid=fopen(logfile);
    tline=fgets(logfid);
    while ischar(tline)
        C=textscan(tline, '%*[^>] %*c %s %*[^\n]');
        if ~isempty(C{1,1})
            ngz_entry=strcat(C{1,1},ext);
            niigz_list=[niigz_list; ngz_entry];
        end
        tline=fgets(logfid);
    end
    fclose(logfid);
end
img_auto_translate(niigz_list,'t',0,0);
%disp(niigz_list);