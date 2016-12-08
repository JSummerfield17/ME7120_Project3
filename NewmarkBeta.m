function [ pos, vel, acc ] = NewmarkBeta( K, M, Beta, gamma, Zeta, deltat, t, tapp )
%%Finds displacement, velocity, acceleration using the newark beta method

clc
close all

run dof_strip
%Reduce stiffness and mass from 3D to 2D by stripping displacement in z, 
%rotation in x and rotation in y. Apply BCS (restrain x and y 
%at node 1 and restrain y at node 51)
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

%Calculate Damping and from Force Matrix
C = Damps(K, Zeta);
R = zeros(size(K,1),1);
R(149,1) = 100000;

time = (0:deltat:t);
T = length(time);
R0 = zeros(size(R,1),size(R,2));

dn = zeros(size(K,1),1); %Initial position zero
ddn = zeros(size(K,1),1); %Initial velocity zero
dd2n = M\R; %F=ma ==  a=F/m

for i = 1:T
    pos(:,i) = dn(:);
    vel(:,i) = ddn(:);
    acc(:,i) = dd2n(:);
    step = time(i);
    %If still within pluck use R, else use R0
    if step <= tapp
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

theta = pos(121,:);
dtheta = vel(121,:);
ddtheta = acc(121,:);
hold on
figure(1)
plot(time,theta) 
ylabel('thetaz41 (rad/s)')
xlabel('time(s)')
figure(2)
plot(time,dtheta) 
ylabel('dthetaz41 (rad/s)')
xlabel('time(s)')
figure(3)
plot(time,ddtheta)
ylabel('ddthetaz41 (rad/s)')
xlabel('time(s)')

%Working Input
%[pos, vel, acc] = NewmarkBeta(K, M, 0.25, 0.5, 0.0001, 0.15, 0.01)

end

