function [ pos, vel, acc ] = NewmarkBeta( K, M, Beta, gamma, deltat, t )
%%Finds displacement, velocity, acceleration using the newark beta method

K = full(K);
M = full(M);
C = Damp(K);

time = [0:deltat:t];
T = size(time,2);
R0 = zeros(size(K,1),1);
R = R0;
R(15,1) = 100000

dn = zeros(size(K,1),1)
ddn = zeros(size(K,1),1)
dd2n = M\R


for i = 1:(T-1)
    step = time(i)
    dn = dn(i)
    ddn = ddn(i)
    dd2n = dd2n(i)
    if step <= t
        dn(i+1) = dnplus1(R, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K);
    else
        dn(i+1) = dnplus1(R0, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K);
    end
    ddn(i+1) = ddnplus1(Beta,deltat,dnplus1,dn,ddn,dd2n,gamma);
    dd2n(i+1) = dd2nplus1(Beta, deltat, dn, ddn, dd2n, dnplus1); 
    pos(i,:,:) = dn(:,:);
    vel(i,:,:) = ddn(:,:);
    acc(i,:,:) = dd2n(:,:);
end

pos
vel
acc
end

