function arrplusone=arrcat(basearr,newstring,MSL)

if nargin<3 
    if isempty(basearr)
        MSL=500;
    else
        max(size(basearr,2),size(newstring,2));
    end
end

arrplusone=[basearr; newstring blanks(MSL-length(newstring))];
return;
end