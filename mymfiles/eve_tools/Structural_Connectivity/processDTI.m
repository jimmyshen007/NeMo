function [ADC,FA,VectorF,EigenValMtx,DifT]=processDTI(DTIdata,H, parameters)
% This function will perform DTI calculations on a given DTI dataset.
% 
% [ADC,FA,VectorF,EigenValMtx]=processDTI(DTIdata,gradients, parameters)
%
% DTIdata:  A 4D matrix containing no-gradient volume (DTIdata(:,:,:,1))
% and gradient volumes (DTIdata(:,:,:,2:end))
% H is the gradients table (ndirs x 3 matrix)
% parameters: A struct containing, parameters.Bvalue (usually 800 or 1000 
%           depending on the image acquisition protocol), parameters.BackgroundTreshold 
%           (suggested value is 400, but this depends on the image values),
%           and parameters.WhiteMatterExtractionThreshold (suggested value 
%           is 0.1), parameters.textdisplay
%           
% ADC: A 3D matrix with the Apparent Diffuse Coefficient (ADC)
% FA: A 3D matrix with the fractional anistropy (FA)
% VectorF: A 4D matrix with the (main) fiber direction in each pixel
% EigenValMtx: A 4D matrix with all eigenvalues for every voxel
% DifT: A 3D matrix with all  Diffusion tensors [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
%
% This function is written by A. Raj Weill-Cornell (June 2009), modified by
% Amy in Dec 2010 to output the Diffusion tensor volume

N1 = [0     1     0;
     1     1     1;
     0     1     0];
N2= [0     0     0;
     0     1     0;
     0     0     0];
nhood = cat(3, N2, N1, N2);

szD = size(DTIdata);
ndirs = szD(4);
if(parameters.textdisplay), disp('Start DTI function'); pause(0.1); end
    
% Find the zero gradient voxel volume(s)
S0=DTIdata(:,:,:,1);
% Make a vector to store the different B-values (timing)
Bvalue=parameters.Bvalue;

% Read the input data (DTIdata), and seperate the zero gradient
% and other gradients into different matrices.
if(parameters.textdisplay), disp('Separate gradient and none gradient datasets'); pause(0.1); end
if(parameters.textdisplay), disp('Create the b matrices'); pause(0.1); end
% Create the b matrices
% (http://www.meteoreservice.com/PDFs/Mattiello97.pdf)
b=zeros([3 3 size(H,1)]);
for i=1:size(H,1),
    b(:,:,i)=Bvalue*H(i,:)'*H(i,:);
end

if(parameters.textdisplay), disp('Voxel intensity to absorption conversion'); pause(0.1); end
% Convert measurement intentensity into absorption (log)
for i=1:ndirs,
    DTIdata(:,:,:,i)=log((DTIdata(:,:,:,i)./(S0+eps))+eps);
end

if(parameters.textdisplay), disp('Create B matrix vector [Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz]'); pause(0.1); end
% Sort all b matrices in to a vector Bv=[Bxx,2*Bxy,2*Bxz,Byy,2*Byz,Bzz];
Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';

% Create a matrix to store the Diffusion tensor of every voxel
% [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
DifT=zeros([size(S0) 6],'single');
% Create a matrix to store the fractional anistropy (FA)
FA=zeros(size(S0),'single');
% Create a matrix to store the Apparent Diffuse Coefficient (ADC)
ADC=zeros(size(S0),'single');

if nargout>2
    % Create a maxtrix to store the (main) fiber direction in each pixel
    VectorF=zeros([size(S0) 3],'single');
end
if nargout>3
    % Create a matrix to store the eigen values of every voxel
    EigenValMtx=zeros([size(S0) 3],'single');
end
if(parameters.textdisplay), disp('Calculate Diffusion tensor, eigenvalues, and other parameters of each voxel'); pause(0.1); end
% Loop through all voxel coordinates
for x=1:size(S0,1),
    for y=1:size(S0,2),
        for z=1:size(S0,3),
            
            % Only process a pixel if it isn't background
            if(S0(x,y,z)>parameters.BackgroundTreshold)
                
                % Calculate the Diffusion tensor [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz],
                % with a simple matrix inverse.
                Z=-squeeze(DTIdata(x,y,z,:));
                %size(Bv(:,2:end)), size(Z(2:end)),
                M=(Bv(2:end,:))\(Z(2:end));
               
                % The DiffusionTensor (Remember it is a symetric matrix,
                % thus for instance Dxy == Dyx)
                DiffusionTensor=[M(1) M(2) M(3); M(2) M(4) M(5); M(3) M(5) M(6)];

                % Calculate the eigenvalues and vectors, and sort the 
                % eigenvalues from small to large
                [EigenVectors,D]=eig(DiffusionTensor); EigenValues=diag(D);
                [t,index]=sort(EigenValues); 
                EigenValues=EigenValues(index); EigenVectors=EigenVectors(:,index);
                EigenValues_old=EigenValues;
                
                % Regulating of the eigen values (negative eigenvalues are
                % due to noise and other non-idealities of MRI)
                if((EigenValues(1)<0)&&(EigenValues(2)<0)&&(EigenValues(3)<0)), EigenValues=abs(EigenValues);end
                if(EigenValues(1)<=0), EigenValues(1)=eps; end
                if(EigenValues(2)<=0), EigenValues(2)=eps; end
                
                % Apparent Diffuse Coefficient
                ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
                
                % Fractional Anistropy (2 different definitions exist)
                % First FA definition:
                %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                % Second FA definition:
                FAv=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                
                % Store the results of this pixel in the volume matrices
                ADC(x,y,z)=ADCv;
                DifT(x,y,z,:)=[DiffusionTensor(1:3) DiffusionTensor(5:6) DiffusionTensor(9)];
                if nargout>3
                    EigenValMtx(x,y,z,:)=EigenValues;
                end
                % Only store the FA and fiber vector of a voxel, if it exceed an anistropy treshold
                if(FAv>parameters.WhiteMatterExtractionThreshold)
                    FA(x,y,z)=FAv;
                    if nargout>2
                        VectorF(x,y,z,:)=EigenVectors(:,end)*EigenValues_old(end);
                    end
                end
            end
        end
    end
end

%% new section to morphologically filter data and create brain mask
L = int8(ADC>0);
L = imerode(L, nhood);
if ~length(find(L))==0
    s  = regionprops(bwlabeln(L), 'Area');
    mxArea = max([s(:).Area]);
    npixremove = round(0.3*mxArea);
    L = bwareaopen(L,npixremove);
    L = imdilate(L, nhood);
end
ADC(L==0) = 0;
FA(L==0) = 0;

if(parameters.textdisplay), disp('DTI function Finished'); pause(0.1); end





