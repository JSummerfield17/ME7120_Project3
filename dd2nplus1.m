function out = dd2nplus1( Rextnplus1, M, Beta, deltat, dn, ddn, dd2n, C, gamma, Keff )
%Find the Dd2n_N+1 result

a = 1/(Beta*deltat^2)*(dnplus1-dn-deltat*ddn);
b = (1/(2*Beta)-1)*dd2n;

dd2nplus1 = a-b;

end

