function result = aKTHunbiased(L,M,upper,psi0)
    cc=load([fileparts(mfilename('fullpath')),'\data\cn2016']);
    dcc=load([fileparts(mfilename('fullpath')),'\data\dcn2016']);
    delta=eye(L);
    Ekr=E(upper,psi0);
    Dn=zeros(upper,1);
    for i=2:upper
        cn=cc.cn2016(i);
        dcn=dcc.dcn2016(i);
        if i>=2&&i<=M
            Dn(i,1)=sigma(i)+dcn;
        elseif i>M
            Dn(i,1)=sigma(i)+cn;
        end
    end
    a=zeros(L,L);
    for i=2:L
        for j=2:L
                EE=Ekr(:,j).*Ekr(:,i);
                EEC=EE(2:upper).*Dn(2:upper);
                sum0=sum(EEC(:));
                a(i,j)=delta(i,j)*Dn(i)-Ekr(j,i)*Dn(j)-Ekr(i,j)*Dn(i)+sum0;
        end
    end
    result=a(2:L,2:L);
end

function result = E(n,psi0)
   e=coe_e(n,psi0);
   Enk=zeros(n,n);
   for i=1:n
       for j=1:n
           Enk(i,j)=e(i,j)*(2*j+1)/2;
       end
   end
   result=Enk;
end

function result = coe_e(n,psi0)
%%------------set parameters-----------------------------------------------
    t = cos(psi0 * pi / 180);
    R = zeros(n+2,n+2);
    P_0 = 1;
    R(1,1) = 1 / 3 * (t^3 + 1);
    P = Pn(n,psi0);
%%-------------calculate---------------------------------------------------
    for i = 2 : 1 : n + 1
        for j = 2 : 1 : n + 1
            if (i ~= j)
                R(i,j)=(i*(i+1)*P(j)*(P(i+1)-P(i-1))/(2*i+1)-j*(j+1)*P(i)*(P(j+1)-P(j-1))/(2*j+1))/((i-j)*(i+j+1));
            end
        end
    end
    R_20=P_0*(P(3)-P(1))/5;
    R(3,1)=(3*(3+1)*P(1)*(P(3+1)-P(3-1))/(2*3+1)-1*(1+1)*P(3)*(P(1+1)-P_0)/(2*1+1))/((3-1)*(3+1+1));
    R(2,2)=9 * R(3,1) / 10 - R_20 / 2 + 3 * R(1,1) / 5;
    for i = 3 : 1 : n + 1 
            R(i,i)=(i+1)*(2*i-1)*R(i+1,i-1)/(i*(2*i+1))-(i-1)*R(i,i-2)/i+(2*i-1)*R(i-1,i-1)/(2*i+1);
    end  
    result = R; 
end

function result = Pn(n,psi0)
    P = zeros(1,n+2);
    t = cos(psi0 * pi / 180);
    P(1) = t;
    P(2) = 3 * t^2 / 2 - 1 / 2;
    for ii = 3 : 1 : n + 2
        P(ii) = (2 * ii - 1 ) * t * P(ii-1) / ii - (ii-1) * P(ii-2) / ii;
    end
    result = P;
end

function result = sigma(n)
    u=0.99899012912;
    C0=1;
    cT=2*C0/u^2;
    sigma=cT*(1-u)*u^n;
    result=sigma;
end