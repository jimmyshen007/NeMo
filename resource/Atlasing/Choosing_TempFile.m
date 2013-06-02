
%=====================Used to be Internal function==================================%
function [TempFile] = Choosing_TempFile(Type)
if nargin==0
    Type = input('Choosing template(T1, T2, PET, SPECT, EPI, PD or Other):   ','s');
    Type = lower(Type);
end
switch Type
    case 't1'
        if exist(fullfile(spm('Dir'),'templates','T1.nii'))
            TempFile = fullfile(spm('Dir'),'templates','T1.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'T1.mnc');
        end
    case 't2'
        if exist(fullfile(spm('Dir'),'templates','T2.nii'))
            TempFile = fullfile(spm('Dir'),'templates','T2.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'T2.mnc');
        end
    case 'pet'
        if exist(fullfile(spm('Dir'),'templates','PET.nii'))
            TempFile = fullfile(spm('Dir'),'templates','PET.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'PET.mnc');
        end
    case 'spect'
        if exist(fullfile(spm('Dir'),'templates','SPECT.nii'))
            TempFile = fullfile(spm('Dir'),'templates','SPECT.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'SPECT.mnc');
        end
    case 'epi'
        if exist(fullfile(spm('Dir'),'templates','EPI.nii'))
            TempFile = fullfile(spm('Dir'),'templates','EPI.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'EPI.mnc');
        end
    case 'pd'
        if exist(fullfile(spm('Dir'),'templates','PD.nii'))
            TempFile = fullfile(spm('Dir'),'templates','PD.nii');
        else
            TempFile = fullfile(spm('Dir') , 'templates' , 'PD.mnc');
        end
    case 'other'
        [TempFileName,TempFilePath] = uigetfile('*.mat;*.img','Reading Reference Template File ...');
        TempFile = [TempFilePath TempFileName];
    otherwise
        disp('Unknown template file ...');
end;
return;

