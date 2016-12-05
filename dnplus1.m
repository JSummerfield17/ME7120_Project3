function dnplus1 = dnplus1(F, M, Beta, deltat, dn, ddn, dd2n, C, gamma, K )
%Find the D_N+1 result
a=((1/(Beta*deltat^2))*M+(gamma/(Beta*deltat))*C+K)^-1;
b=a*F;
c=a*(M*((1/(Beta*deltat^2))*dn+(1/(Beta*deltat))*ddn+(1/(2*Beta)-1)*dd2n));
d=a*(C*((gamma/(Beta*deltat))*dn+(gamma/Beta-1)*ddn+(gamma/Beta-2)*(deltat/2)*dd2n));
dnplus1=(b+c+d);

end


