function out = ddnplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, C, gamma, Keff )
%Find the Ddn_N+1 result

a = (gamma/(Beta*deltat))*(dnplus1-dn);
b = (gamma/beta-1)*ddn;
c = deltat*(gamma/(2*Beta)-1)*d2dn;

ddnplus1 = a-b-c;

end

