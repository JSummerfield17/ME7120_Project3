function [C] = C_local( K_local, M_local, Zeta )
%Finds local damping

C = Zeta*(2*M_local*sqrt(K_local/M_local));


end

