function mpc = matsysLSHg(Wp,mpc,sol,xinit)
% matsysLSHg sets up the stacked matrices L and S as well as Hessian H and 
% gradient g

%% Get wind turbine system and input matrices, MPC parameters
mpc = matsys(Wp,mpc);

L = mpc.AA;
S = mpc.BBt;
Nh = mpc.Nh;
N = Wp.turbine.N;

%% Caluclate S_tilde, L_tilde, 

% Calculate S_tilde and L_tilde
Isum   = kron(eye(Nh),ones(1,N));
Csel   = kron(eye(Nh*N),[0 1 0]);
C_tilde = Isum * Csel;
L_tilde = C_tilde * L;
S_tilde = C_tilde * S;

% Set utemp
ukneg1 = sol.turbine.CT_prime;
utemp = repmat(ukneg1,Nh,1); % + randn(Nh*N,1); 

% Calculate gradient of dU'dU: 2*[dU1-dU2, dU2-dU3,...dUN]
% or alternatively dU = deltaU* U * + dC
deltaUvec = utemp - [ukneg1; utemp(1:(N*(mpc.Nh-1)))];
dDeltaUdU = eye(N*Nh) + [zeros(N,N*Nh);- eye(N*(Nh-1)), zeros(N*(Nh-1),N)];
dCU = [-eye(N); zeros(N*Nh-N,N)];
deltaUvec0 =  dCU * ukneg1;

% deltaUvec1 = dDeltaUdU * utemp + dCU * ukneg1; %for tests

% Calculate error trajectory (with and without shift by S_tilde* utemp)
E0 = (L_tilde*xinit - mpc.Pref(sol.k:sol.k+Nh-1));
E = E0 + S_tilde* utemp; % (L_tilde*xinit + S_tilde* utemp - mpc.Pref(sol.k:sol.k+Nh-1));

%% Compute Hessian H, derivative g = dJk/dUk %%
% min Jk(Uk) s.t. h1(Uk) = 0, h2(Uk) < 0, at current estimate xl
% <=> min( x'H(xl)x + g(xl)x s.t. dh1(xl)x+h1(xl)=0, dh2(xl)x+h2(xl)< 0
% with 
% Jk = E'Q E+(Uk-Uss)'R(Uk-Uss) + Psi
% dJk/dUk = 2S'Q (Lx+SUk-Xss) + 2R(Uk-Uss)

% Gradient with shift: mpc.g = 2*(S_tilde'*mpc.Q*E + mpc.R*grad_dU)';     
mpc.g = 2*(E' * mpc.Q * S_tilde + deltaUvec' * mpc.R * dDeltaUdU);

% Gradient without shift: mpc.g = 2*(S_tilde'*mpc.Q*(Lx0 -Pref) + mpc.R*grad_dU)';
mpc.g0 = 2*(E0' * mpc.Q * S_tilde + deltaUvec0' * mpc.R * dDeltaUdU);

% Hessian: H = 2*(S_tilde'*mpc.Q*S_tilde + mpc.R*RH); 
H = 2*(S_tilde' * mpc.Q * S_tilde + dDeltaUdU' * mpc.R * dDeltaUdU); 
mpc.H = 1/2*(H + H');

%% Assign outputs
mpc.L_tilde = L_tilde;
mpc.S_tilde = S_tilde;
mpc.utemp = utemp;

