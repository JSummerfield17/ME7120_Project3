function C = Damp( K )
%Finds Damping matrix from Zeta, K, M

Zeta = .02; 

%Determine Mode Shapes and Eigenvalues
[mode, lam ]= eig(K);
w = sqrt(lam);

lam_mat = 2*Zeta*w;

C = (mode')^-1*lam_mat*mode^-1;
end

