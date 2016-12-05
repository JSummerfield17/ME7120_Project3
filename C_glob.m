function [Cg] = C_glob( K, M, Zeta )
%Builds global damping matrix

A = size(K,1);
B = size(K,2);
Cg = zeros(A,B);

for i = 1:A
    for j = 1:B
        Cg(i,j) = C_local(K(i,j), M(i,j), Zeta);
    end
end

end

