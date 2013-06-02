function [Lmax Nmin] = calcula_Lmax(DirSignal)
%%%%%%%
pos = ismember(DirSignal, [0 0 0], 'rows');
ind = find(pos);
DirSignal(ind,:) = [];

N = length(DirSignal);
Lmax = 0;
while (Lmax+1)*(Lmax+2)/2 <= N
  Lmax = Lmax +2;
end
Lmax = Lmax - 2;
Nmin = sum(1:4:2*Lmax+1);
