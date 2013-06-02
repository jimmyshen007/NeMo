function [V,F] = getFV(maxlevel,shape,r)

if nargin == 1
    shape = 'ico';
    r = 1;
end
if nargin == 2
    r = 1;
end    
FV = sphere_tri(shape,maxlevel,r);

V = FV.vertices;
F = FV.faces;