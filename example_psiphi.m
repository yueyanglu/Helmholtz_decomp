

%% create a vel fld

clear

% coord
[JDM,IDM] = deal(100,100); % 3min for [100,1000]
[xemin, xemax] = deal(-10, 10);           
[yemin, yemax] = deal(-10, 10); 
[xc,xe,dx] = FDGrid(xemin,xemax,JDM); 
[yc,ye,dy] = FDGrid(yemin,yemax,IDM);

% C-grid mesh
[plon,plat] = meshgrid(xc,yc);         % p-grid [JDM   - IDM  ]
[qlon,qlat] = meshgrid(xe,ye);         % q-grid [JDM+1 - IDM+1]
[ulon,ulat] = meshgrid(xe,yc);         % u-grid [JDM   - IDM-1]
[ulon,ulat] = deal(ulon(:,2:end-1),ulat(:,2:end-1));
[vlon,vlat] = meshgrid(xc,ye);         % v-grid [JDM-1 - IDM  ]
[vlon,vlat] = deal(vlon(2:end-1,:),vlat(2:end-1,:));

% size of the grid
[njp, nip] = size(plon);
[njq, niq] = size(qlon);
[nju, niu] = size(ulon);
[njv, niv] = size(vlon);

% mesh size of C-grid
[cpx,cpy] = deal(1/dx * ones(njp,nip), 1/dy * ones(njp,nip)); 
[cqx,cqy] = deal(1/dx * ones(njq,niq), 1/dy * ones(njq,niq)); 
[cux,cuy] = deal(1/dx * ones(nju,niu), 1/dy * ones(nju,niu)); 
[cvx,cvy] = deal(1/dx * ones(njv,niv), 1/dy * ones(njv,niv)); 

cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);


%% create the velocity fld to be decomposed

% VEL_psi (rotational comp) (-y/r, x/r)
u_psi = - ulat ./ sqrt(ulon.^2 + ulat.^2);
v_psi =   vlon ./ sqrt(vlon.^2 + vlat.^2);

% VEL_phi (divergent comp) (x/r^2, y/r^2)
u_phi = ulon ./ sqrt(ulon.^2 + ulat.^2);
v_phi = vlat ./ sqrt(vlon.^2 + vlat.^2);

[u,  v] = deal(u_psi + u_phi,  v_psi + v_phi);

%% do

alpha = 1e-16;

whichmethod = 2; % 1 for MATLAB's fminunc; 2 for minFunc

%---------------------------------------------- options for optmz

if whichmethod == 1
    opt = optimoptions('fminunc');
    opt.Algorithm = 'quasi-newton'; % trust-region, quasi-newton
    opt.HessUpdate = 'bfgs';
    opt.SpecifyObjectiveGradient = true;
    opt.Display = 'iter';
    opt.MaxFunEvals = 1e6;
    opt.MaxIter = 1e3;
    opt.OptimalityTolerance = 1e-10;
    opt.DerivativeCheck = 'off'; % check if supplied grad match FD approximations
    opt.FiniteDifferenceType = 'central';
else
    opt.Method = 'lbfgs';
    opt.GradObj = 'on';
    opt.Display = 'iter';
    opt.MaxFunEvals = 1e6;
    opt.MaxIter = 1e6;
    opt.optTol = 1e-10; % Tolerance on the first-order optimality
    opt.DerivativeCheck = 'off';
    opt.numDiff = 0;
end

%---------------------------------------------- initial guess of psi/phi
psi0 = randn(njq,niq);
phi0 = randn(njp,nip);

%---------------------------------------------- psi/phi that minimize ja
tic;
[psi_ot,phi_ot] = uv2psiphi(psi0,phi0,u,v,cxy,alpha,opt,whichmethod);
toc;

%---------------------------------------------- psi/phi to uv

dpsidy = (psi_ot(2:end,2:end-1) - psi_ot(1:end-1,2:end-1)) .* cuy; % u-
dpsidx = (psi_ot(2:end-1,2:end) - psi_ot(2:end-1,1:end-1)) .* cvx; % v-

dphidy = (phi_ot(2:end,:) - phi_ot(1:end-1,:)) .* cvy; % v-
dphidx = (phi_ot(:,2:end) - phi_ot(:,1:end-1)) .* cux; % u-

% components of flux
[u_psi_re, v_psi_re] = deal( - dpsidy, dpsidx);
[u_phi_re, v_phi_re] = deal(   dphidx, dphidy);
[u_re, v_re] = deal(u_psi_re + u_phi_re, v_psi_re + v_phi_re);



