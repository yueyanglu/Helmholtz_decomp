function adj = derive_adj(err,cxy,sz,ZBC,MBC,idNoN)
% 
% Compute the adjoint term in the gradient of Tikhionov's functional, ie.,
% the first term in eqn (11).
% Formulation:  adj = A' * (y - A*x)
% 
% This function essentially calculates the curl and div of the velocity
% difference, ie., the error matrix.
% 
% adj = A' * vel = [- du/dy + dv/dx] = [ curl of vel ]
%                  [  du/dx + dv/dy]   [ div of vel  ]
% 
% The adj term has the same size with psi/phi vec [njp*nip + njq*niq,  1].
% 

% size of the grid
[njp, nip] = deal(sz.njp, sz.nip);
[njq, niq] = deal(sz.njq, sz.niq);
[nju, niu] = deal(sz.nju, sz.niu);
[njv, niv] = deal(sz.njv, sz.niv);

% mesh factors when calc finite difference
cpx = cxy.cpx;
cpy = cxy.cpy;
cqx = cxy.cqx;
cqy = cxy.cqy;

% reshape error mat, same size with full y (flx vector with NaN)
errful = zeros(nju*niu + njv*niv, 1);
errful(idNoN) = err;

% 2d flx difference, on u-/v- points
u = reshape(errful(1 : nju*niu), [nju niu]);
v = reshape(errful(nju*niu+1 : nju*niu+njv*niv), [njv niv]);


%------------------------------------------------------ curl and div of flx
% curl and div on q-/p- point, respectively. Note to take care of BC
curl = 0 * zeros(njq,niq);
div = 0 * zeros(njp,nip);

% differentiation terms of uflx and vflx. But smaller size than psi and phi
dvdx = (v(:,2:end) - v(:,1:end-1)) .* cqx(2:end-1,2:end-1); % njq-2, niq-2
dudy = (u(2:end,:) - u(1:end-1,:)) .* cqy(2:end-1,2:end-1); % njq-2, niq-2
% 
dudx = (u(2:end-1,2:end) - u(2:end-1,1:end-1)) .* cpx(2:end-1,2:end-1); % njp-2, nip-2
dvdy = (v(2:end,2:end-1) - v(1:end-1,2:end-1)) .* cpy(2:end-1,2:end-1); % njp-2, nip-2

% curl (q-) & div (p-)
curl(2:end-1,2:end-1) = dvdx - dudy;
div(2:end-1,2:end-1) = dudx + dvdy;

% BC?
if strcmpi(ZBC,'closed') && strcmpi(MBC,'closed')
    
end
    
% reshape to form adj term, same size with x
adj = [reshape(curl,[njq*niq 1]); reshape(div,[njp*nip 1])];



