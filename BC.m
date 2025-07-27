function X = BC(X,lb,ub,dim)
for j=1:dim
    if X(j)<lb(j)
        X(j) = lb(j)+rand*(ub(j)-lb(j));
    end
    if X(j)>ub(j)
        X(j) = lb(j)+rand*(ub(j)-lb(j));
    end   
end  
end