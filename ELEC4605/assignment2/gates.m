%% const is either alpha beta
%% array : either xj or km
function [U_gates,Dum] = gates(n,dx, const,l,k,const2)
      L = mod(l-1,n);  % L and K are index of bits |LK>, l and k are counters for LK
      K = mod(k-1,n);
      fprintf("%d%d\n",L,K);
      dum_var = 1; % counter for Dum
      for jl = 1:2
         for jk = 1:2
%              U_PE(jl,jk,idx) = expm(-1j*beta*(dx^2)*(mod(jl,2)*2^(l)+lambda)*(mod(jk,2)*2^(k)+lambda));
               Dum(dum_var) = (mod(jl+1,2)*2^(l)+const2)*(mod(jk+1,2)*2^(k)+const2);
               dum_var = dum_var + 1;
               fprintf("jl djk = %d%d\n",mod(jl+1,2),mod(jk+1,2));
         end
      end
      U_gates = zeros(4,4);
      for i = 1:4
          U_gates(i,i) = exp(-1j*const*(dx.^2)*Dum(i));
      end
      fprintf("preparation for 4x4 matrix completed!\n");
end