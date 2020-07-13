function [psi,phi] = uv2psiphi(psi0,phi0,uflx,vflx,cxy,alpha,opt,whichmethod)
%UV2PSIPHI computes the streamfunction and potential from given vectors.
% 
%   [psi,phi] = UV2PSIPHI(psi0,phi0,uflx,vflx,cxy,alpha,opt) performs the 
%   computation by solving a minimization problem proposed by Li et al
%   (2006). 
%   The method has advantages in extracting the rotational (psi-related)
%   and divergent (phi-related) components of the vector field without 
%   specifying boundary conditions and computational efficiencies.
% 
%   The first dimension of variables are latitude (y-direction) for
%   convenience. Typical Arakawa C grid is used except that the W/E edges 
%   of u and S/N edges of v points are deleted. So niu = nip - 1 and
%   njv = njp -1. This is consistent with the function 'derive_A.m'. 
%    
%   Input:
%      psi0 [njq,niq]: Initial guess of streamfunction
%      phi0 [njp,nip]: Initial guess of potential
%      uflx [nju,niu]: Zonal flx 
%      vflx [njv,niv]: Meridional flx
%      cxy           : Struct that saves the factors for finite differece
%      alpha         : Tikhonov's regularization parameter
%      opt           : Optimization options
%
%   Output:
%      psi  [njp,nip]: Streamfunction
%      phi  [njp,nip]: Potential
% 
% 
%   Note 1: U and V surrounding land should be masked with NaNs. Initial
%           psi and phi should be masked with 0.
% 
%   Note 2: flx = flx_psi + flx_phi
%               = k X del_psi + del_phi
%           Thus,
%           uflx = - dPsi/dy + dPhi/dx
%           vflx =   dPsi/dx + dPhi/dy
%     
%  Author: Yueyang Lu, 2020.
% 
%  Sources: 
%  Li, Z., Y. Chao, and J. C. McWilliams, 2006: Computation of the
%  streamfunction and velocity potential for limited and irregular domains.
%  Mon. Wea. Rev., 134, 3384–3394
%  The python code based on the same method can be found in (Tiago Bilo)
%  https://github.com/iuryt/vector_fields. Note that the configure of grid 
%  and basic logic of the code are very different.
% 

%%

%-------------------------------------------------------- size of vars
[nju,niu] = size(uflx);
[njv,niv] = size(vflx);
[njq,niq] = size(psi0);
[njp,nip] = size(phi0);
% put in a structure
sz = struct('njp',njp,'nip',nip,'njq',njq,'niq',niq,'nju',nju,'niu',niu,...
    'njv',njv,'niv',niv);

%---------------------------------------------------------- col vectors

% initial psi and phi vec 'x'
x0 = [reshape(psi0,[njq*niq 1]); reshape(phi0,[njp*nip 1])];

% flx vec 'y'
y = [reshape(uflx,[nju*niu 1]); reshape(vflx,[njv*niv 1])];
% set NaNs to 0
isnanY = isnan(y);
y(isnanY) = 0;


%------------------------------------------------------------ Minimization

%------ prepare params for optimization

% 1. matrix A to be sent to func, since remain constant in iteration
A = derive_A(cxy,sz);

% Set rows of A to zero according to Nans in u/v, i.e., enforces 0 = 0 in 
% the linear system.
A(isnanY,:) = 0;

% 2. handle of the obj func & its gradient (Jacobian)
fh = @(x)objfun(x,y,A,alpha);


%------ perform the optmz
if whichmethod == 1
    fprintf('Optimization using MATLAB''s ''fminunc'' ...\n')
    [x_optm, fval] = fminunc(fh,x0,opt);
else
    fprintf('Optimization using ''minFunc'' ...\n')
    [x_optm, fval] = minFunc(fh,x0,opt);
end
fprintf('F(x): %f\n',fval)


%------ reshape to get 2d psi & phi
psi = reshape(x_optm(1 : njq*niq), [njq niq]);
phi = reshape(x_optm(njq*niq+1 : njq*niq+njp*nip), [njp nip]);


%% objective function & its Jacobian

function [f,gradf] = objfun(x,y,A,alpha)
% Obj func and gradient. Eqn (10) and (11) in Li et al (2006),respectively.
        
    % err vector
    err = y - A * x;
    
    % Tikhionov's functional [1-by-1].
    f = .5 * (err' * err) + .5 * alpha * (x' * x);
    
    % Gradient (Jacobian) of the functional at each qp-point [#qp-by-1].
    gradf = - A' * err + alpha .* x; 
    
end

end
