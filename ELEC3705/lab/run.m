%% Initialisation
down = [1;0];
up = [0;1];
psi = (up+down)/sqrt(2);
psi = up;
%% Run the following commands to generate a mesh and draw the Husimi Q distribution
mesh = gen.sphereMesh();
coherent = qTop.coherentStateSphere(mesh,1/2);
Q = qTop.analysis.Husimi(psi, coherent);
qTop.analysis.plotHusimi(mesh, Q);
