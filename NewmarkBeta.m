function [ pos, vel, acc ] = NewmarkBeta( K, M, F, Beta, gamma, deltat, t )
%%Finds displacement, velocity, acceleration using the newark beta method

clc

K = full(K);
K = K(1:30,1:30);
M = full(M);
M = M(1:30,1:30);
C = Damp(K);

time = (0:deltat:t);
T = size(time,2);
R = full(F(1:30,1));
R0 = zeros(size(R,1),size(R,2));

dn = zeros(size(K,1),1); %Initial position zero
ddn = zeros(size(K,1),1); %Initial velocity zero
size(M);
size(R);
dd2n = M\R; %F=ma


for i = 1:T
    pos(i,:) = dn(:);
    vel(i,:) = ddn(:);
    acc(i,:) = dd2n(:);
    step = time(i);
    %If still within pluck use R, else use R0
    if step <= t
        Dn = dnplus1(R, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K);
    else
        Dn = dnplus1(R0, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K);
    end
    Ddn = ddnplus1(Beta,deltat,Dn,dn,ddn,dd2n,gamma);
    Dd2n = dd2nplus1(Beta, deltat, dn, ddn, dd2n, Dn);
    dn = real(Dn);
    ddn = real(Ddn);
    dd2n = real(Dd2n);
end

end

