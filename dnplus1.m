function dnplus1 = dnplus1(Rnplus1, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K )
%Find the D_N+1 result
a=(1/(Beta*deltat^2))*M+(gamma/(Beta*deltat))*C+K;
b=Rnplus1;
c=(M*((1/(Beta*deltat^2))*dn+(1/(Beta*deltat))*ddn+(1/(2*Beta)-1)*dd2n));
d=(C*((gamma/(Beta*deltat))*dn+(gamma/Beta-1)*ddn+(gamma/Beta-2)*(deltat/2)*dd2n));
dnplus1t =a\(b+c+d);
dnplus1 = real(dnplus1t);
end


