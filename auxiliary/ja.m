function J_alpha = ja(x,y,cxy,sz,idNoN,alpha)
% 
% Calculate the scalar Tikhionov's functional, ie, the sum of the objective 
% functional and Tikhonov's regularization term, as equation (10) in 
% Li et al (2006).
% 
% J_alpha(x) = J(x) + J_reg(x)
%            = .5 * (y - Ax)' * (y - Ax)  +  .5 * alpha * x' * x
% 

% Derive flx from PSI and PHI, i.e., tht operation: A * x
Ax = derive_Ax(x,cxy,sz,idNoN);

% Error matrix 
err = y - Ax;

% objective functional
J = .5 * (err' * err);

% Tikhonov's regularization term
J_reg = .5 * alpha * (x' * x);

% Tikhionov's functional
J_alpha = J + J_reg;
