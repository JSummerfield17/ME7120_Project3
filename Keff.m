function Kf = Keff( Beta, deltat, M, gamma, C, K )
%Calculates K_effective to then be used in the calculation of D_N+1

a = 1/(Beta*deltat^2)*M;
b = gamma/(Beta*deltat)*C;

Kf = a+b+K;

end

