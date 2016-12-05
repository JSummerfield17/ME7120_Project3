function [ pos, vel, acc ] = Newmark_Beta( K, C, M, Beta, gamma, deltat, t )
%Finds displacement, velocity, acceleration using the newark beta method

time = linspace(0,t,deltat);
T = size(time,1);

dn = zeros(size(K,1),1);
ddn = zeros(size(K,1),1);
dd2n = zeros(size(K,1),1);


t = .01;
deltat = 0.0001;
time = [0:deltat:t];
T = length(time);
F = zeros(T,1);

for i = 1:T
    if time(i)>=0 & time(i)<=t;
        F(i) = 100000;
    else F = 0;
    end
end

pos = zeros(T,3);

for i = 1:(T-1)
    dn = dn(i);
    ddn = ddn(i);
    dd2n = dd2n(i);
    dn(i+1) = dnplus1(F, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K);
    ddn(i+1) = ddnplus1(Beta,deltat,dnplus1,dn,ddn,dd2n,gamma);
    dd2n(i+1) = dd2nplus1(Beta, deltat, dn, ddn, dd2n, dnplus1); 
    pos(i,:,:) = dn(:,:);
    vel(i,:,:) = ddn(:,:);
    acc(i,:,:) = dd2n(:,:);
end


end

