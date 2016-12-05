%Inputs
%K, M 
Zeta=0.02;%damping
t1=0;%initial time
deltat=0.01;%step size
t2=0.1;%end time
impactime=0.01;%impact time
F=100000;
Beta=1/4;%Beta
gamma=1/2;%gamma
initialdis=0;
initialvel=0;

E=1;
I=1;
A=1;
L=1;
G=1;
J=1;
f=(6*L);
g=(4*L*L);
h=(2*L*L);
k=(E*A)/(L);
m=(G*J)/(L);
%element 1 is the beam element in xy plane
 K =[12 f -12 f;
     f g -f h;
     -12 -f 12 -f;
      f h -f g];
  P=100;
M=zeros(4,4);
M(1,1)=2;
M(2,2)=2;
M(3,3)=2;
M(4,4)=2;
M(1,3)=1;
M(2,4)=1;
M(:,1)=M(1,:);
M(:,2)=M(2,:);
M(:,3)=M(3,:);
M(:,4)=M(4,:);
M=((P*A*L)/6)*M



Cg=zeros(4,4);
dt=deltat;

% Cg = C_glob( K, M, Zeta );
% Kf = Keff( Beta, deltat, M, gamma, Cg, K );
% dn=zeros(size(K,1),1);%initial disp assuming to be zero
% ddn=zeros(size(Cg,1),1);%initial vel assuming to be zero
% %dn(1,1)=0;%
% %ddn(1,1)=0;%
% dd2n=F*M^-1; %F=ma
for T=t1:dt:t2;
dn=zeros(size(K,1),1);
ddn=zeros(size(Cg,1),1);
dn(1,1)=initialdis;
ddn(1,1)=initialvel;
Rextnplus1=zeros(size(K,1),1);
Rextnplus1(1,1)=F;
% Calculate Dppn0 using F=ma;
dd2n=Rextnplus1*M^-1
dd2n=dd2n'
Rextnplus1=zeros(size(K,1),1) 
    
if T==impactime
    Rextnplus1=zeros(size(K,1),1);
   F(T)=F;
   Rextnplus1(1,1)=F(T);
    dnplus1 = dnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, K );
    ddnplus1 = ddnplus1(Beta,deltat,dnplus1,dn,ddn,dd2n,gamma);
    dd2nplus1 = dd2nplus1(Beta, deltat, dn, ddn, dd2n, dnplus1);

else 
   
    Rextnplus1=zeros(size(K,1),1);
    dnplus1 = dnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, Cg, gamma, K );
    ddnplus1 = ddnplus1(Beta,deltat,dnplus1,dn,ddn,dd2n,gamma);
    dd2nplus1 = dd2nplus1(Beta, deltat, dn, ddn, dd2n, dnplus1);

end
end
plot(dnplus1,T)
    
   