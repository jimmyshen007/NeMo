function [ output_args ] = StringOut(StringIn,fid )
%LineOutToFile: simply uses fprintf without having to worry about fid's and other silly syntax

fprintf(fid,'%s\n',StringIn);

end

