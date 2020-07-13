function A = derive_A(cxy,sz)
% 
% Derive the matrix A that calc u/v from psi and phi by finite difference. 
% The result of A*x using A is EXACTLY the same with 'derive_Ax.m'.
% A is made up of 4 block matrices having different physical meanings.
% 
% The C-grid used is the typical one excaept that the W/E edges of u and
% S/N edges of v are deleted. This exempts us from specifying BCs.
% 
% cxy  -- struct that saves factors for finite difference
% szs  -- struct that saves the size of vars on the C-grid
% 

% size of the grid
[njp, nip] = deal(sz.njp, sz.nip);
[njq, niq] = deal(sz.njq, sz.niq);
[nju, niu] = deal(sz.nju, sz.niu);
[njv, niv] = deal(sz.njv, sz.niv);

% mesh factors when calc finite difference
% cpx = cxy.cpx;
% cpy = cxy.cpy;
% cqx = cxy.cqx;
% cqy = cxy.cqy;

cux = cxy.cux;
cuy = cxy.cuy;
cvx = cxy.cvx;
cvy = cxy.cvy;

%%  
%---------------------------------------------- A11, get u from PSI, -d/dy
[J,I,L1d] = deal([]);
mm = 0; % linear index of non-zeros entries in Aji
    
% loop over psi-points
for i = 1:niq
    for j = 1:njq
        
        % linear index of psi-vec, also the col-id for Aji
        m = sub2ind([njq,niq], j,i);
                
        %------------ 1. Look to the SOUTH, standing on psi(j,i)
        % subscripts of u-points on the S of psi(j,i)
        [ju,iu] = deal(j-1,i-1);
        
        if ju >= 1 && iu >= 1 && ju <= nju && iu <= niu
            % 
            mm = mm + 1;
            
            % linear index of u-vec, also the row-id for Aji
            n = sub2ind([nju,niu], ju,iu);
            
            % J for u, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = - cuy(ju,iu);
        end
        
        %------------ 2. Look to the NORTH
        % subscripts of u-points on the N of psi(j,i)
        [ju,iu] = deal(j,i-1);
        
        if ju >= 1 && iu >= 1 && ju <= nju && iu <= niu
            %
            mm = mm + 1;
            
            % linear index of u-vec, also the row-id for Aji
            n = sub2ind([nju,niu], ju,iu);
            
            % J for u, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = cuy(ju,iu);
        end
        
    end
end

% A11 is a sparse matrix [nju*niu -by- njq*niq].
A11 = sparse(J,I,L1d,nju*niu,njq*niq);


%---------------------------------------------- A21, get v from psi, d/dx
[J,I,L1d] = deal([]);
mm = 0; % linear index of non-zeros entries in Aji
    
% loop over psi-points
for i = 1:niq
    for j = 1:njq
        
        % linear index of psi-vec, also the col-id for Aji
        m = sub2ind([njq,niq], j,i);
                
        %------------ 1. Look to the WEST
        % subscripts of v-points on the W of psi(j,i)
        [jv,iv] = deal(j-1,i-1);
        
        if jv >= 1 && iv >= 1 && jv <= njv && iv <= niv
            % 
            mm = mm + 1;
            
            % linear index of v-vec, also the row-id for Aji
            n = sub2ind([njv,niv], jv,iv);
            
            % J for v, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = cvx(jv,iv);
        end
        
        %------------ 2. Look to the EAST
        % subscripts of v-points on the E of psi(j,i)
        [jv,iv] = deal(j-1,i);
        
        if jv >= 1 && iv >= 1 && jv <= njv && iv <= niv
            %
            mm = mm + 1;
            
            % linear index of v-vec, also the row-id for Aji
            n = sub2ind([njv,niv], jv,iv);
            
            % J for v, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = - cvx(jv,iv);
        end
        
    end
end

% A21 is a sparse matrix [njv*niv -by- njq*niq].
A21 = sparse(J,I,L1d,njv*niv,njq*niq); 


%------------------------------------------------ A12, get u from PHI, d/dx
[J,I,L1d] = deal([]);
mm = 0; % linear index of non-zeros entries in Aji
    
% loop over phi-points
for i = 1:nip 
    for j = 1:njp
        
        % linear index of phi-vec, also the col-id for Aji
        m = sub2ind([njp,nip], j,i);
                
        %------------ 1. Look to the WEST, standing on phi(j,i)
        % subscripts of u-points on the W of phi(j,i)
        [ju,iu] = deal(j,i-1);
        
        if iu >= 1
            % 
            mm = mm + 1;
            
            % linear index of u-vec, also the row-id for Aji
            n = sub2ind([nju,niu], ju,iu);
            
            % J for u, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = cux(ju,iu);
        end
        
        %------------ 2. Look to the EAST
        % subscripts of u-points on the E of phi(j,i)
        [ju,iu] = deal(j,i);
        
        if iu <= niu
            %
            mm = mm + 1;
            
            % linear index of u-vec, also the row-id for Aji
            n = sub2ind([nju,niu], ju,iu);
            
            % J for u, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = - cux(ju,iu);
        end
        
    end
end

% A12 is a sparse matrix [nju*niu -by- njp*nip].
A12 = sparse(J,I,L1d,nju*niu,njp*nip); 


%------------------------------------------------ A22, get v from PHI, d/dy
[J,I,L1d] = deal([]);
mm = 0; % linear index of non-zeros entries in Aji

% loop over phi-points
for i = 1:nip
    for j = 1:njp
        
        % linear index of phi-vec, also the col-id for Aji
        m = sub2ind([njp,nip], j,i);
        
        %------------ Look to the SOUTH
        % subscripts of v-points on the S of phi(j,i)
        [jv,iv] = deal(j-1,i);
        
        if jv >= 1
            %
            mm = mm + 1;
            
            % linear index of v-vec, also the row-id for Aji
            n = sub2ind([njv,niv], jv,iv);
            
            % J for v, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = cvy(jv,iv);
        end
        
        %------------ Look to the NORTH
        % subscripts of v-points on the N of phi(j,i)
        [jv,iv] = deal(j,i);
        
        if jv <= njv
            % 
            mm = mm + 1;
            
            % linear index of v-vec, also the row-id for Aji
            n = sub2ind([njv,niv], jv,iv);
            
            % J for v, I for phi
            [J(mm), I(mm)] = deal(n,m);
            L1d(mm) = - cvy(jv,iv);
        end
        
    end
end

% A22 is a sparse matrix [njv*niv -by- njp*nip].
A22 = sparse(J,I,L1d,njv*niv,njp*nip); 

% form A [#_uv -by- #_qp]
A = [ A11 A12
    A21 A22 ];


