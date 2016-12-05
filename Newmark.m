%Inputs
%K, M 
Zeta=0.02;%damping
F=[100000 0 0 0];%force just for trial run
t1=0;%initial time
deltat=0.01;%step size
t2=0.1;%end time
impactime=0.01;%impact time
Beta=1/4;%Beta
gamma=1/2;%gamma

Cg = C_glob( K, M, Zeta );
Kf = Keff( Beta, deltat, M, gamma, Cg, K );
dn=zeros(size(K,1),1);%initial disp assuming to be zero
ddn=zeros(size(Cg,1),1);%initial vel assuming to be zero
%dn(1,1)=0;%
%ddn(1,1)=0;%
dd2n=F*M^-1; %F=ma
for T=t1:deltat:t2       
if T==impactime
    Rextnplus1=F;
    dnplus1 = dnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );
    ddnplus1 = ddnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );
    dd2nplus1 = dd2nplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );

else 
   
    Rextnplus1=0;
    dnplus1 = dnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );
    ddnplus1 = ddnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );
    dd2nplus1 = dd2nplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, Kf );

end
end
plot(dnplus1,T)

    
   