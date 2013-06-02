function trk_amass(range,qbdir,isofn,new_trk_out,atlassize,cleanup)

%trk_amass.m

trk_list={};

for m=range(1):range(2)
    trk_list{m}=qbelong_to_trk(num2str(m),qbdir,isofn,atlassize);
end

trk_list=char(trk_list);

if nargin<6
    concat_trk(trk_list,new_trk_out);
else
    concat_trk(trk_list,new_trk_out,cleanup);
end
end