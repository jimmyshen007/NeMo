% [Ylm] = Ylm(l,m,theta,phi)
%
% Spherical Surface Harmonic 
%
% These spherical coordinate conventions are the same as those used by
% the Mathematica SphericalHarmonicY[l,m,theta,phi] command.

% $Id: Ylm.m,v 1.4 2005/09/05 12:10:00 rodgers Exp $

% Copyright: Chris Rodgers, Aug, 2005.
% Modified : Erick Canales-Rodriguez, 2005.
%
% This file is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This file is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this file; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% See the accompanying LICENSE file for details.

function Ylm=Ylm(l,m,theta,phi)

if l<0 || abs(m)>l || fix(m)~=m || fix(l)~=l || isempty(theta) || isempty(phi)
  error('Spherical order out of range!')
  return;
end

% Calculate the harmonic
tmp = legendre(l,cos(reshape(theta,1,[])));

absm=abs(m);
if m>=0
  Ylm = reshape(tmp(m+1,:),size(theta)).*exp(i*m*phi);
else
  Ylm = ((-1)^absm * factorial(l-absm)/factorial(l+absm)) * reshape(tmp(absm+1,:),size(theta)).*exp(i*m*phi);
end

Ylm = sqrt((2*l+1)*factorial(l-m)/(4*pi*factorial(l+m)))*Ylm;
