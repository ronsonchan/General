% 
function [control_0, control_1] = decomposition(n,U_gates)
% decompose 4x4 into 2 2x2.
control_0= eye(2);
control_1= eye(2); 
i = 1;
k = 1;
for j = 1:4
   if(mod(j+1,2) == 0)
       control_0(i,i) = U_gates(j,j);
       i = i+1;
   else
       control_1(k,k) = U_gates(j,j);
       k = k+1;
   end
end
fprintf("2 level decomposition completed!\n");
end