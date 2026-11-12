function mpc = HSqLd_WF(sol,Wp,mpc,ll,xinit) %N,q,R_,p,xk,Xk,DT,ref,us)
% HSqLd calculate matrices L, S, H and g for Lemke's algorithm
%
% Inputs: * N: number of samples optimization (N = 10)
%         * q: diag Q state weighting matrix stage cost l(x,u) = x'Qx + u'Ru
%         * R_: R_ = kron(eye(N),R), R input weighting matrix stage cost
%         * xk: quasiLPV KF state estimates
%         * Xk: Xko(= L*xk+S*U) delayed by one sample timestep
%         * Ts: sample time
%         * ref: Reference signal
%         * us: steady state input
%
% Outputs: * L: stacked matrix L for state trajectory Xk = L x_k + S U_k
%          * S: stacked matrix S for state trajectory
%          * H: H essian: H = 2*(S'*Q*S+R_)
%          * g: Derivative dJk/dUk: g = 2*S'*Q*(L*xk-ref) - 2*(R_*Usp);

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WF simulation: Antje Dittmer

nx = 3; %length(xk);  
N =  Wp.turbine.N;
ni = N; %;
Nh = mpc.Nh;
DT = Wp.sim.h;

mpc.cf = sol.turbine.cf;
mpc.cp = sol.turbine.cp;

if ll==1
    rhoP = repmat(sol.turbine.Ur,1,mpc.Nh);
else
    rhoP = mpc.V;
end
    
Phi = sol.turbInput.phi *pi/180;

%%%%%%%%% Build L, S matrices Xk = L*xk+S*Uk%%
[Aj,Bj] = ABqLd_WF(DT,rhoP(:,1),Phi,mpc);

nxWF = nx*N; %number of states windfarm;
L = zeros(nxWF*Nh,nxWF);
S = zeros(nx*N,ni*N);
L(1:nxWF,1:nxWF) = Aj;
S(1:nxWF,1:ni) = Bj;

for idx = 1:Nh-1
    [Aj,Bj] = ABqLd_WF(DT,rhoP(:,idx+1),Phi,mpc);  
    L(idx*nxWF+1:(idx+1)*nxWF,1:nxWF) = Aj*L((idx-1)*nxWF+1:idx*nxWF,1:nxWF);
    S(idx*nxWF+1:(idx+1)*nxWF,1:(idx+1)*ni) = [Aj*S((idx-1)*nxWF+1:idx*nxWF,1:idx*ni),Bj];
end

%% Calculate S_tilde and L_tilde
Isum   = kron(eye(Nh),ones(1,N));
Csel   = kron(eye(Nh*N),[0 1 0]);
C_tilde = Isum*Csel;
L_tilde = C_tilde*L;
S_tilde = C_tilde*S;

% Calculate utemp
if ~isfield(mpc, 'U')
    X = L*xinit + S*repmat(sol.turbine.CT_prime,Nh,1);
    Y = X;
    utemp = Y(mpc.MU);
else
    utemp = mpc.U;
end

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
deltaUvec = utemp - [sol.turbine.CT_prime; utemp(1:(N*(mpc.Nh-1)))];
grad_dU = deltaUvec - [deltaUvec(N+1:end);zeros(N,1)];

% calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
R2 = [zeros(N*(Nh-1),N),- eye(N*(Nh-1)); zeros(N,N*Nh)];
RH = diag([2*ones(N*(Nh-1),1);ones(N,1)])+ R2 + R2';
E1 = (L_tilde*xinit + S_tilde* utemp - mpc.Pref(sol.k:sol.k+Nh-1));

% Gradient:-
mpc.g = 2*(S_tilde'*mpc.Q*E1 + mpc.R*grad_dU)'; % max(abs(g(:)- mpc.g(:))./g(:))

% Hessian:-
H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH);
mpc.H = 1/2*(H + H'); % max(abs(H(:)- mpc.H(:))./H(:))

% temp = mu_old
mpc.utemp = utemp;
mpc.AA  = L;
mpc.BBt = S;
mpc.CC  = kron(eye(Nh),blkdiag(mpc.c{:}));

end

