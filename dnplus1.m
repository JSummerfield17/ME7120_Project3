function [ output_args ] = dnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, C, gamma, Keff )
%Find the D_N+1 result
a = Rextnplus1;
b = M*[(1/(Beta*deltat^2)*dn)+1/(Beta*deltat)*ddn+(1/(2*Beta)-1)*d2dn];
c = C*[(gamma/(Beta*deltat^2)*dn+(gamma/Beta-1)*ddn+deltat*(gamma/(2*Beta)-1)*d2dn];

d = a+b+c;

ddnplus1 = d\Keff;

end
end

