function [psi,u_psi,v_psi,phi,u_phi,v_phi,output] = uv_decomp(u,v,cxy)
% uv [JDM,IDM]
% 

%% check the size of inputs

if ~isequal(numel(u),numel(v),numel(cxy.cux),numel(cxy.cvx))
    error('Size of the inputs must be the same !!!');
end

[JDM,IDM] = size(u);

%% prepare the shape of vars to satisfy 'uv2psiphi.m'

% CHOOSE subregion, with the W/E edges of u and S/N edges of v are deleted
[jjq, iiq] = deal(1:JDM, 1:IDM);                % E.g. 6-6
[jjp, iip] = deal(jjq(1:end-1), iiq(1:end-1));  % E.g. 5-5
[jju, iiu] = deal(jjq(1:end-1), iiq(2:end-1));  % E.g. 5-4
[jjv, iiv] = deal(jjq(2:end-1), iiq(1:end-1));  % E.g. 4-5

% vars
u = u(jju,iiu);
v = v(jjv,iiv);
[cqx,cqy] = deal(cxy.cqx(jjq,iiq), cxy.cqy(jjq,iiq)); 
[cpx,cpy] = deal(cxy.cpx(jjp,iip), cxy.cpy(jjp,iip)); 
[cux,cuy] = deal(cxy.cux(jju,iiu), cxy.cuy(jju,iiu)); 
[cvx,cvy] = deal(cxy.cvx(jjv,iiv), cxy.cvy(jjv,iiv)); 

% size of the grid
[njq, niq] = deal(length(jjq),length(iiq));
[njp, nip] = deal(length(jjp),length(iip));
[nju, niu] = deal(length(jju),length(iiu));
[njv, niv] = deal(length(jjv),length(iiv));

cxy = struct('cpx',cpx,'cpy',cpy,'cqx',cqx,'cqy',cqy,'cux',cux,'cuy',cuy,...
    'cvx',cvx,'cvy',cvy);

% disp 
fprintf('Input u/v size %s\n', mat2str([JDM IDM]))
fprintf('Q-point size %s\n', mat2str([njq niq]))
fprintf('P-point size %s\n', mat2str([njp nip]))
fprintf('U-point size %s\n', mat2str([nju niu]))
fprintf('V-point size %s\n', mat2str([njv niv]))

%% do

% Tikhionov's regularization parameter
alpha = 1e-14;

% Which optimz function to use. 1 for MATLAB's fminunc; 2 for minFunc
whichmethod = 2; 

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
    opt.Display = 'final';
    opt.MaxFunEvals = 1e8;
    opt.MaxIter = 1e8;
    opt.optTol = 1e-12; % Tolerance on the first-order optimality
    opt.DerivativeCheck = 'off';
    opt.numDiff = 0;
end

%---------------------------------------------- initial guess of psi/phi
psi0 = rand(njq,niq);
phi0 = rand(njp,nip);
disp('Random initial guess of psi and phi !!')

%---------------------------------------------- psi/phi that minimize ja
[psi_ot,phi_ot,output] = uv2psiphi(psi0,phi0,u,v,cxy,alpha,opt,whichmethod);


%% put results onto original grid

[psi,u_psi,v_psi,phi,u_phi,v_phi] = deal(NaN * zeros(JDM,IDM));

%---------------------------------- psi and phi
psi(jjq,iiq) = psi_ot;
phi(jjp,iip) = phi_ot;

%---------------------------------- components of u, v
dpsidy = (psi_ot(2:end,2:end-1) - psi_ot(1:end-1,2:end-1)) .* cuy; % u-
dpsidx = (psi_ot(2:end-1,2:end) - psi_ot(2:end-1,1:end-1)) .* cvx; % v-

dphidy = (phi_ot(2:end,:) - phi_ot(1:end-1,:)) .* cvy; % v-
dphidx = (phi_ot(:,2:end) - phi_ot(:,1:end-1)) .* cux; % u-

% optimal rot and div comp
[u_psi(jju,iiu), v_psi(jjv,iiv)] = deal( - dpsidy, dpsidx);
[u_phi(jju,iiu), v_phi(jjv,iiv)] = deal(   dphidx, dphidy);


