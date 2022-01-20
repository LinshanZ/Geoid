function result = bRCRbiased(L,M,n,sigma0)
    QL=QLn(L,n,sigma0);
    bn=zeros(M,1);
    for i=2:M
        bn(i)=-QL(i);
    end
    result=bn(2:M);
end