function mpc = HSqLdD_WF(sol,Wp,mpc,ll,xinit,mpctest) %N,q,R_,p,xk,Xk,DT,ref,us)
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

nx = mpc.nx;
nxV = nx*2; % For now as all state are included in the output
N = Wp.turbine.N;
ni = N; %; one input per turbine
Nh = mpc.Nh; % history
DT = Wp.sim.h;
nxWF = nx * N;
nxV_WF = nxV*N; %number of states windfarm;

mpc.cf = sol.turbine.cf;
mpc.cp = sol.turbine.cp;

if ll==1
    rhoP = repmat(sol.turbine.Ur,1,mpc.Nh);
else
    rhoP = mpc.V;
end

Phi = sol.turbInput.phi *pi/180;

%%%%%%%%% Build L, S matrices Xk = L*xk + S*Uk%%
[AV,BV] = ABqLdD_WF(DT,rhoP(:,1),Phi,mpc);
%[AVa,BVa] = ABqLd_WF(DT,rhoP(:,1),Phi,mpc,1);
% Ad = AV(nxWF+1:end,nxWF+1:end); % sys matrix original system
% Bd = BV(nxWF+1:end,:); % input matrix original system
Ad2 = AV(1:nxWF,nxWF+1:end);
%AdV = AV; % sys matrix original system
%BdV = BV; % input matrix original system

LambdaV = zeros(nxV_WF*Nh,nxV_WF);
SV = zeros(nxV*N,ni*N);
LambdaV(1:nxV_WF,1:nxV_WF) = AV;
SV(1:nxV_WF,1:ni) = BV;

for idx = 1:Nh-1
    [AV,BV] = ABqLdD_WF(DT,rhoP(:,idx+1),Phi,mpc);
    LambdaV(idx*nxV_WF+1:(idx+1)*nxV_WF,1:nxV_WF) = AV*LambdaV((idx-1)*nxV_WF+1:idx*nxV_WF,1:nxV_WF);
    SV(idx*nxV_WF+1:(idx+1)*nxV_WF,1:(idx+1)*ni) = [AV*SV( (idx-1)*nxV_WF+1 : idx*nxV_WF, 1:idx*ni),BV];
end

[Ad1,Bd1] = ABqLd_WF(DT,rhoP(:,1),Phi,mpctest);
Lambda = zeros(nxWF*Nh,nxWF);
S = zeros(nxWF*Nh,ni*Nh);
Lambda(1:nxWF,1:nxWF) = Ad1;
S(1:nxWF,1:ni) = Bd1;

for idx = 1:Nh-1
    [Ad1,Bd1] = ABqLd_WF(DT,rhoP(:,idx+1),Phi,mpc); % For LPV
    Lambda(idx*nxWF+1:(idx+1)*nxWF,1:nxWF) = ...
        Ad1*Lambda((idx-1)*nxWF+1:idx*nxWF,1:nxWF);
    S(idx*nxWF+1:(idx+1)*nxWF,1:(idx+1)*ni) = ...
        [Ad1*S((idx-1)*nxWF+1:idx*nxWF,1:idx*ni),Bd1];
end

%% Calculate S_tilde and L_tilde

Isum = kron(eye(Nh),ones(1,N));
Csel1 = kron(eye(N),[0 1 0]); % select turbine power states at sample k
% Csel = kron(eye(Nh),Csel1); % select turbine power states at all samples
% C_tilde = Isum*Csel;
% L_tilde = C_tilde*Lambda;
% S_tilde = C_tilde*S;

% Calculate utemp
% if ~isfield(mpc, 'U')
%X = Lambda*xinit + S*repmat(sol.turbine.CT_prime,Nh,1);
utemp = repmat(sol.turbine.CT_prime,Nh,1); % X(mpc.MU);
% else
%     utemp = mpc.U;
% end

% States [F,P,cthead,deltaF,deltaP, deltaCthead]
Csel1V   = [Csel1, zeros(size(Csel1))]; % select turbine power states at sample k
CselV   = kron(eye(Nh),Csel1V); % select turbine power states at all samples
CV_tilde = Isum*CselV;
L_D = CV_tilde*LambdaV;
S_D = CV_tilde*SV;

% Calculate utemp
%mpc_MUV = repmat([mpc.MU(1:nxWF); zeros(nxWF,1)],Nh,1);
%xinitD = [xinit; zeros(3*N,1)];

Ad = mpctest.AA(1:nxWF,1:nxWF);
Bd = mpctest.BBt(1:nxWF,1:ni);
temp = [xinit; Ad2\((Ad-eye(size(Ad)))*xinit + Bd*sol.turbine.CT_prime)];
mpc.xinitD = temp;

% if ~isfield(mpc, 'U')
%     X1 = LambdaV*mpc.xinitD;
%     utemp1 = X1(mpc_MUV == 1);
% else
%     utemp1 = mpc.U;
% end
% %utemp = utemp1;

% test
% x1Test = mpctest.AA(1:nxWF,1:nxWF) * xinit + mpctest.BBt(1:nxWF,1:ni)* sol.turbine.CT_prime; % Ad * xinit + Bd * sol.turbine.CT_prime;
% x2Test = mpctest.AA(1:nxWF,1:nxWF) * x1Test + mpctest.BBt(1:nxWF,1:ni)* sol.turbine.CT_prime; % Ad * xinit + Bd * sol.turbine.CT_prime;
% %x2aTest = mpctest.AA(1:nxWF,1:nxWF) ^2 * xinit + [mpctest.AA(1:nxWF,1:nxWF)*mpctest.BBt(1:nxWF,1:ni) + mpctest.BBt(1:nxWF,1:ni) ]* sol.turbine.CT_prime; % Ad * xinit + Bd * sol.turbine.CT_prime;
% x2aTest = mpctest.AA(1:nxWF,1:nxWF) ^2 * xinit + [mpctest.AA(1:nxWF,1:nxWF)*mpctest.BBt(1:nxWF,1:ni) , mpctest.BBt((1:nxWF)+ nxWF,(1:ni)+ni) ]* utemp(1:ni*2); % Ad * xinit + Bd * sol.turbine.CT_prime;
% y2bTest = L_tilde(2,:) * xinit + S_tilde(2,1:ni*2)* utemp(1:ni*2);
% 
% Csel1 * [mpctest.AA(1:nxWF,1:nxWF)*mpctest.BBt(1:nxWF,1:ni) , mpctest.BBt((1:nxWF)+nxWF,(1:ni)+ni)]
% 
% yTest = Csel1 * x1Test;
% yTest2 = Csel1 * x2Test;
% yTest2a = Csel1 * x2aTest;
% sum(yTest2a)

% x1TestV = AdV * mpc.xinitD;
% 
% x2TestV = LambdaV((1:nxV_WF)+nxV_WF,:) * mpc.xinitD;
% yTestV = Csel1V * x1TestV;
% yTestV2 = Csel1V * x2TestV;
% 
% sum(yTestV2)

%% Gradient and Hessian: DeltaU:velocity form input and gradient dDeltaU/dU
deltaUvec = utemp - [sol.turbine.CT_prime; utemp(1:(N*(Nh-1)))];
dDeltaUdU = eye(N*Nh) + [zeros(N,N*Nh);- eye(N*(Nh-1)), zeros(N*(Nh-1),N)];


% xTraj = Lambda * xinit  + S * utemp;
% xTrajD = LambdaV * mpc.xinitD + SV * deltaUvec;
% 
% Ytraj = mpctest.L_tilde * xinit  + mpctest.S_tilde * utemp;
% Ytraj_V = L_D * mpc.xinitD;
% Error trajectory of velocity form
% E =  (mpctest.L_tilde * xinit  + mpctest.S_tilde * utemp - mpc.Pref(sol.k:sol.k+Nh-1));
ED = (L_D * mpc.xinitD  + S_D * deltaUvec - mpc.Pref(sol.k:sol.k+Nh-1));

%% Gradient and Hessian
mpc.g = 2*(ED'* mpc.Q * S_D * dDeltaUdU + deltaUvec' * mpc.R * dDeltaUdU);
H = 2*(dDeltaUdU' * S_D' * mpc.Q * S_D * dDeltaUdU ...
    + dDeltaUdU' * mpc.R * dDeltaUdU);
mpc.H = (H+H')/2;


% temp = mu_old
mpc.utemp = utemp;
mpc.AA  = LambdaV;
mpc.BBt = SV;
mpc_cV = [eye(nxWF),zeros(nxWF)];
mpc.CC  = kron(eye(Nh),mpc_cV); %kron(eye(Nh),blkdiag(mpc.c{:}));

mpc.L_D = L_D;
mpc.S_D = S_D;

return

mpc.Ulim_lower = 0.2 *ones(Nh*N,1); %#ok<UNRCH>
mpc.Ulim_upper = 2*ones(Nh*N,1);
Ac0_ = [- eye(Nh*N); eye(Nh*N)];
bc0_ = [- mpc.Ulim_lower; mpc.Ulim_upper];
bc = bc0_ - Ac0_ * mpc.utemp;

%     Ac0_1 = [-tril(ones(Nh*N)); tril(ones(Nh*N))];
%     bc0_1 = [- mpc.Ulim_lower + sol.turbine.CT_prime ;
%         mpc.Ulim_upper - sol.turbine.CT_prime];
%% QP solution for standard form
% Calculate the solution to quadprog uTilde. This code
%%
% # Calculates the solution of the quadratic sub-problem $\tilde{U}$
% # Calculates the optimal control input as U = $\tilde{U}$+ $U^l$ and based
% on it the cost function $J(U)$
% # Calculated the corresponding $\Delta U$ and checks that the cost function
% $J_V(\Delta U) = J(U)$

Aeq = []; beq = [];lb = [];ub = [];x0 = mpc.utemp;
opts = optimoptions('quadprog'); opts.Display = 'none';
[uTilde,fval,EXITFLAG] = quadprog(mpctest.H, mpctest.g, Ac0_,bc,Aeq,beq,lb,ub,x0,opts);
costQPorig = 0.5* uTilde' *mpc.H* uTilde + mpc.g*uTilde; %Check fval

uOpt = uTilde + mpc.utemp; % utilde = uOpt - utemp: Get ctrl input of orig problem
X = L_tilde * xinit + S_tilde * uOpt; % state trajectory based on uopt
Etemp_L = mpc.Pref(sol.k:sol.k+Nh-1) - X; % tracking error value trajectorz
d_utemp2 = [uOpt(1:N) - sol.turbine.CT_prime;...
    uOpt(N+1:end) - uOpt(1:end-N)]; % change in input value
% Calculate the cost function
mpctest.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;

deltaUvecOpt2 = uOpt - [sol.turbine.CT_prime; uOpt(1:(N*(Nh-1)))]; %delta U from U
XD2 = L_D * mpc.xinitD + S_D * deltaUvecOpt2;
Etemp_LD2    = mpc.Pref(sol.k:sol.k+Nh-1) - XD2; % tracking error value
% Cost function
costValue2 = Etemp_LD2'*mpc.Q*Etemp_LD2 + deltaUvecOpt2'*mpc.R*deltaUvecOpt2;

%% QP solution for velocity form
% Calculate the solution to quadprog uTilde. This code
%%
% # Calculates the solution of the quadratic sub-problem $\tilde{U}$ fo $J_V(U)$.
% It confirms that the derived solution gives a smaller cost for this QP problem
% than inserting the original solution
% # Calculates the optimal control input as U = $\tilde{U}$+ $U^l$ and based
% on it the cost function $J_V(U)$It confirms that the derived solution gives
% a larger cost for  than inserting the original solution - this is unexpected.
% # Calculates the optimal control input as  $\Delta U = \Delta \tilde{U}$ and
% based on it the cost function $J_V(\Delta U)$

opts.OptimalityTolerance = 1e-8;
[uTilde1,fval1,EXITFLAG1] = quadprog(mpc.H,mpc.g,Ac0_,bc,Aeq,beq,lb,ub,x0,opts);
costQP = 0.5* uTilde1' * mpc.H* uTilde1 + mpc.g*uTilde1;
costQP2 = 0.5* uTilde' * mpc.H* uTilde + mpc.g*uTilde;


end

