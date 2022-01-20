function result = bKTHoptimal(L,M,upper,psi0,k)
    load([fileparts(mfilename('fullpath')),'\data\cn2016'],'cn16');
    cn=cn16;
    load([fileparts(mfilename('fullpath')),'\data\dcn2016'],'dcn16');
    dcn=dcn16;
    b = zeros(M,1);
    sn = sTSVD(L,M,upper,psi0,k);
    QL=QLn(L,M,upper,psi0,k);
    for i = 2 : M
        b(i) = (sn(i-1)+QL(i))*cn(i)/(cn(i)+dcn(i));
    end
    result=b(2:M);
end