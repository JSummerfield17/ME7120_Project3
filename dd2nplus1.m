function dd2nplus1 = dd2nplus1(Beta, deltat, dn, ddn, dd2n, dnplus1)
%Find the Dd2n_N+1 result

a = 1/(Beta*deltat^2)*(dnplus1-dn-deltat*ddn);
b = (1/(2*Beta)-1)*dd2n;

dd2nplus1 = a-b;

end

