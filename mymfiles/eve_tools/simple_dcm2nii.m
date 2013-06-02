function imgout=simple_dcm2nii(indir,inifile,atopt)
startup_varsonly;
%T2src='/home/eve/Documents/T2_DCM';
%dirlist=dir(T2src);
if nargin < 3
    atopt=0;
    if nargin<2
        inifile=[start_dir filesep 'eve_tools' filesep 'T1_dcm2nii.ini'];
    end
end
logfile=[start_dir filesep 'd2n_log.txt'];

imgout=[];


    %indir=[T2src filesep dirlist(i).name];
disp(['system call to ''dcm2nii'': Converting ' indir ' ...']);
system([start_dir filesep 'mricron/dcm2nii -b ' inifile ' ' indir ' > ' logfile]);
    
logfid=fopen(logfile);
tline=fgets(logfid);
while ischar(tline)
    C=textscan(tline, '%*[^>] %*c %s %*[^\n]');
    if ~isempty(C{1,1})
        out_entry=strcat('co',C{1,1});
            %out_entry=strcat(C{1,1},'.img');
        imgout=[imgout; out_entry];
        break;
    end
    	tline=fgets(logfid);
end

fclose(logfid);

if atopt
    img_auto_translate(imgout{1,1});
end

disp(imgout{1,1});

return;
end
