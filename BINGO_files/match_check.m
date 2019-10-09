function check=match_check(X,Y)

%This function checks if the sizes of the given cell arrays match, that is,
%there is the same number of cells, and that the second dimensions of the
%cells match (corresponding to time points in the data). Its output is 0,
%if the sizes match, and 1 if they don't.

check=0;
if norm(size(X)-size(Y))>.5
    check=1;
else
    sX=zeros(size(X,2),1);
    sY=zeros(size(X,2),1);
    for j=1:size(X,2)
        sX(j)=size(X{j},2);
        sY(j)=size(Y{j},2);
    end
    if norm(sX-sY)>.1
        check=1;
    end
end




