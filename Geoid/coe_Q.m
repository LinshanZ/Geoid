function result = coe_Q(n,psi0)
%%------------set parameters-----------------------------------------------
    t = cos(psi0 * pi / 180);
    P_0=1;
    U = zeros(1,n+2);
    V = zeros(1,n+2);
    Qn = zeros(1,n);
    upper = 10000;
    P = Pnn(upper,psi0);
    enk = coe_e(upper,psi0);
%%-------------calculate---------------------------------------------------    
    if psi0 ~= 0 
        if t ~=0
            U(1)=log(2/(1-t+sqrt(2-2*t)));
            V0=log(1+2/sqrt(2-2*t));
        else
            U(1)=0;
            V0=0;
        end
        V(1)=t*V0+sqrt(2-2*t)-1;
        U(2)=t*U(1)-sqrt(2-2*t)-P(1);
        V(2)=(3*t*V(1)-V0+sqrt(2-2*t))/2;
        U(3)=(3*t*U(2)-U(1)-sqrt(2-2*t)+(P_0-P(2))/3)/2;
        for nn=4:1:n+2
            U(nn)=((2*nn-3)*t*U(nn-1)-(nn-2)*U(nn-2)-sqrt(2-2*t)+(P(nn-3)-P(nn-1))/(2*nn-3))/(nn-1);
        end
        for nn=3:1:n+2
            V(nn)=((2*nn-1)*t*V(nn-1)-(nn-1)*V(nn-2)+sqrt(2-2*t))/nn;
        end
        U1 = U(1);
        U2 = U(2) + 1;
        U3 = U(3) + 1 / 2 + t;
        for i = 2 : n
            Pn=P(i);
            Pnjia1=P(i+1);
            Pnjian1=P(i-1);
            Un=U(i)+1/(i-1);
            Unjia1=U(i+1)+1/i+t/(i-1);
            Unjia2=U(i+2)+1/(i+1)+t/i+(3*t^2-1)/(2*(i-1));
            Vnjian1=V(i-1)-1/i-t/(i+1)-(3*t^2-1)/(2*(i+2));
            Vn=V(i)-1/(i+1)-t/(i+2);
            Vnjia1=V(i+1)-1/(i+2);
            Qn(i)=i*(i+1)/((2*i+1)*(i-1)*(i+2))*(Pn*(2*(2*i+1)/(i*(i+1))*(U1-U3)-...
                (i+2)*(Un-Unjia2)-(i-1)*(Vnjia1-Vnjian1))+(Pnjia1-Pnjian1)*(3*U2-(i+2)*Unjia1...
                +(i-1)*Vn))-(2*i^2+2*i+1)/((i-1)*(2*i+1)^2)*Pn*(Pnjia1-Pnjian1)+(2*i+1)/(i-1)*enk(i,i);
        end
    else
        for ii = 2 : 1 : n
            Qn(ii) = 0;
            for jj = 2:1:upper
                f=(2*jj+1)*enk(ii,jj)/(jj-1);
                Qn(ii) = Qn(ii) + f;
            end
        end
    end
    result = Qn;
end

function result = coe_e(n,psi0)
%%------------set parameters-----------------------------------------------
    t = cos(psi0 * pi / 180);
    R = zeros(n+2,n+2);
    P_0 = 1;
    R(1,1) = 1 / 3 * (t^3 + 1);
    P = Pnn(n,psi0);
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

function result = Pnn(n,psi0)
    t = cos(psi0 * pi / 180);
    P = zeros(1,n+2);
    P(1) = t;
    P(2) = 3 * t^2 / 2 - 1 / 2;
    for ii = 3 : 1 : n + 2
        P(ii) = (2 * ii - 1 ) * t * P(ii-1) / ii - (ii-1) * P(ii-2) / ii;
    end
    result = P;
end