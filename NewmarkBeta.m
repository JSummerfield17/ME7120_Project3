function [ pos, vel, acc ] = NewmarkBeta( K, M, Beta, gamma, deltat, t )
%%Finds displacement, velocity, acceleration using the newark beta method

clc

run dof_strip
bcs = [1 2 152];
K = full(K);
K(stripdof,:) = [];
K(:,stripdof) = [];
K(bcs,:) = [];
K(:,bcs) = [];
M = full(M);
M(stripdof,:) = [];
M(:,stripdof) = [];
M(bcs,:) = [];
M(:,bcs) = [];
C = Damp(K);
R = zeros(size(K,1),1);
R(150,1) = 100000;

time = (0:deltat:t);
T = length(time);
R0 = zeros(size(R,1),size(R,2));

dn = zeros(size(K,1),1); %Initial position zero
ddn = zeros(size(K,1),1); %Initial velocity zero
dd2n = M\R; %F=ma ==  a=F/m

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

