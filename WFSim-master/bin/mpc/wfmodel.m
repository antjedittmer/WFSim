function mpc = wfmodel(sol,Wp,mpc,ll)
% wfmodel builds input and system matrices for wind farm model

% build B matrix of the turbine models for current times step
initB = ~isfield(mpc,'bcoef');
uPhi = sol.turbine.Phi * pi/180 * 0;

for kk = 1: Wp.turbine.N
 
    uPhiRep =  repmat(uPhi(kk),mpc.Nh,1)*0;
    if ll==1
        
        mpc.cf = sol.turbine.cf;
        mpc.cp = sol.turbine.cp;
        
        v = repmat(sol.turbine.Ur(kk),mpc.Nh,1);    
    else
        v = mpc.V(kk,:);
    end
    
    if initB
        mpc.bcoef{kk} = mpc.num(2)*[ -mpc.cf(kk); mpc.cp(kk); 1]; % changing init with zeros
    end
    
    % mpc.b{kk} = [mpc.bcoef{kk}(1)*v(1)^2 mpc.bcoef{kk}(2)*v(1)^3 mpc.bcoef{kk}(3)]';
    for nn = 1:mpc.Nh
        mpc.bt(:,kk,nn) = [mpc.bcoef{kk}(1)*(v(nn)*cos(uPhiRep(nn)))^2 mpc.bcoef{kk}(2)*(v(nn)*cos(uPhiRep(nn)))^3 mpc.bcoef{kk}(3)]';
    end
    
end

% build wind farm model
mpc.A = blkdiag(mpc.a{:});
% mpc.B = blkdiag(mpc.b{:});
for kk = 1:mpc.Nh
    for nn = 1:Wp.turbine.N
        b{nn}     = mpc.bt(:,nn,kk);
    end
    mpc.Bt(:,:,kk) = blkdiag(b{:});
end
mpc.C = blkdiag(mpc.c{:});

end

