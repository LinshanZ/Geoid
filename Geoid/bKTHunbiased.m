function result = bKTHunbiased(L,M,upper,psi0,k)
    b = zeros(M,1);
    sn = sTSVD(L,M,upper,psi0,k);
    QL=QLn(L,M,upper,psi0,k);
    for i = 2 : M
        b(i) = sn(i-1)+QL(i);
    end
    result=b(2:M);
end