function ddnplus1 = ddnplus1(Beta,deltat,dnplus1,dn,ddn,dd2n,gamma)
%Find the Ddn_N+1 result

a = (gamma/(Beta*deltat))*(dnplus1-dn);
b = ((gamma/Beta)-1)*ddn;
c = deltat*((gamma/(2*Beta))-1)*dd2n;

ddnplus1t = a-b-c;
ddnplus1 = real(ddnplus1t);

end

