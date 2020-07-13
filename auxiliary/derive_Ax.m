function Ax = derive_Ax(x,cxy,sz,idNoN)
% 
% Derive fluxes from psi and phi fields, ie., y = A*x. No BC is applied
% in this operation!
% 
% flx = flx_psi + flx_phi
%     = k X del_psi + del_phi
% Thus,
% uflx = - dpsidy + dphidx
% vflx =   dpsidx + dphidy
% 
% See equation (2) in Li et al (2006).
% 

% size of the grid
[njp, nip] = deal(sz.njp, sz.nip);
[njq, niq] = deal(sz.njq, sz.niq);
[nju, niu] = deal(sz.nju, sz.niu);
[njv, niv] = deal(sz.njv, sz.niv);

% mesh factors when calc finite difference
cux = cxy.cux;
cuy = cxy.cuy;
cvx = cxy.cvx;
cvy = cxy.cvy;

% reshape to 2d psi & phi
psi = reshape(x(1 : njq*niq), [njq niq]);
phi = reshape(x(njq*niq+1 : njq*niq+njp*nip), [njp nip]);

%-------------------------------------------------- calc flx from psi & phi 
% flx on u-/v- point

% differentiation of streamfunc and potential
dpsidy = (psi(2:end,2:end-1) - psi(1:end-1,2:end-1)) .* cuy; % on u-point
dpsidx = (psi(2:end-1,2:end) - psi(2:end-1,1:end-1)) .* cvx; % v-

dphidy = (phi(2:end,:) - phi(1:end-1,:)) .* cvy; % v-
dphidx = (phi(:,2:end) - phi(:,1:end-1)) .* cux; % u-

% flx
uflx = - dpsidy + dphidx;
vflx =   dpsidx + dphidy;

% Organize uv
Ax = [reshape(uflx,[nju*niu 1]); reshape(vflx,[njv*niv 1])];

% Remove NaNs in flx
Ax = Ax(idNoN);



