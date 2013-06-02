function V=strip_b0(DTIseries)

[d,f,e]=fileparts(deblank(DTIseries));

V1=spm_vol(DTIseries);
b0vol=V1(1);

b0vol.fname=[d filesep f '_b0' e];

V2=spm_read_vols(V1);
b0dat=V2(:,:,:,1);

spm_write_vol(b0vol,b0dat);

V=b0vol.fname;

openMRIcron(b0vol.fname,0);
end