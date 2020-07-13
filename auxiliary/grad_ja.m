function gja = grad_ja(x,y,cxy,sz,idNoN,ZBC,MBC,alpha)
% 
% Calculate the gradient of the Tikhionov's functional, ie., the Jacobian 
% of J_alpha, as equation (11) in Li et al (2006).
% Formulation:  gja = - A' * (y - A*x) + alpha * x
% 
% Ax computes the velocity from the psi and phi, that is, y = A * x is
% essentially 
% [ u ]   =  [ - d/dy   d/dx] * [ psi ]
% [ v ]      [   d/dx   d/dy]   [ phi ].
% 
% So the transpose (adjont matrix) of A is still  
% A' =  [ - d/dy   d/dx] 
%       [   d/dx   d/dy] 
% 
% Thus, the adjont is actually curl and conv:
% - A' * vel = [   du/dy - dv/dx] = [ - curl of vel ]
%              [ - du/dx - dv/dy]   [ - divg of vel ]
%              

% Derive velocity from PSI and PHI, i.e., the operation: A * x
Ax = derive_Ax(x,cxy,sz,idNoN);

% Error vec. Note NaNs in y and Ax have been deleted. 
err = y - Ax;

% Compute the adjoint term, i.e., curl and div of the flux difference 
adj = derive_adj(err,cxy,sz,ZBC,MBC,idNoN);

% Jacobian (gradient) of J_alpha (Tikhionov's functional)
gja = - adj + alpha .* x;