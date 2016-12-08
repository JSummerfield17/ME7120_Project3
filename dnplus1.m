function dnplus1 = dnplus1(Rnplus1, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K )
%Find the D_N+1 result
a=(1/(Beta*deltat^2))*M+(gamma/(Beta*deltat))*C+K
b=Rnplus1;
size(M);
size(dn);
size(ddn);
size(dd2n);
c=(M*((1/(Beta*deltat^2))*dn+(1/(Beta*deltat))*ddn+(1/(2*Beta)-1)*dd2n));
d=(C*((gamma/(Beta*deltat))*dn+(gamma/Beta-1)*ddn+(gamma/Beta-2)*(deltat/2)*dd2n));
size(a);
size(b);
size(c);
size(d);
test = (b+c+d);
size(test);
dnplus1t =a\(b+c+d);
dnplus1 = real(dnplus1t);
end


