function [CT_prime,phi,mpc] = qLMPCcontroller(sol,Wp,mpc)
% qLMPCcontroller calculates optimal MPC CT_prime values based on current
% effective wind speeds and thrust coeffecients via QP

%% Parameter for optimization
Nh = mpc.Nh;        % prediction trajectory length
N = Wp.turbine.N;   % number turbines

% QP boundary condition: (Aco_) x <= bco_
Ac0_ = [- eye(Nh*N); eye(Nh*N)];
bc0_ = [- mpc.Ulim_lower; mpc.Ulim_upper];

% xinit = [Fk Pk CT_prime] 
xinit         = zeros(mpc.nx*Wp.turbine.N,1);
xinit(mpc.Mf) = sol.turbine.force;
xinit(mpc.Mp) = sol.turbine.power;
xinit(mpc.Mu) = sol.turbine.CT_prime;

% nl-1 is number of times the rotor-averaged wind speeds in the horizon
% will be updated during one sample. If nl=1, the rotor-averaged wind speeds
% are taken constant in the horizon
nl = 2;
Uopt0 = NaN(Nh*N,nl);

Aeq = [];
beq = [];
lb = [];
ub = [];
x0 = [];
options = optimoptions('quadprog','Display','none');

for ll = 1:nl
    
    % build wind farm input and system matrices
    mpc  = wfmodel(sol,Wp,mpc,ll);
    
    % build Tpelitz matrices for prediction horizon, Hessian, and gradient 
    mpc  = matsysLSHg(Wp,mpc,sol,xinit);
    
    % solve QP optimization problem and calculate time
    tic
    bc   = bc0_ - Ac0_ * mpc.utemp; % A(u_tilde + utemp) <= Bco
    [uval,fval] = quadprog(mpc.H,mpc.g,Ac0_,bc,Aeq,beq,lb,ub,x0,options);
    mpc.x = uval + mpc.utemp;
    mpc.uval = uval;
    mpc.fval = fval;
    mpc.displaytime = toc;
    
    %% True wind state space model of the wind farm
    uTilde = mpc.uval + mpc.utemp;
    X    = mpc.AA*xinit + mpc.BBt*uTilde;
    Y    = mpc.CC*X;                                %yopt
    P    = reshape(Y(mpc.MP),N,Nh);  %Power output
 
    % tracking error value
    Etemp_L    = mpc.Pref(sol.k:sol.k+Nh-1) - sum(P)';
    
    % change in input error value
    d_utemp2 = [uTilde(1:N) - sol.turbine.CT_prime;...
        uTilde(N+1:end) - uTilde(1:end-N)];
    
    % Calculate the cost function
    mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
    
    Uopt0(:,ll) = uTilde; %value(Y(mpc.MU))
    Popt       = P;
    
    temp       = (repmat(Popt(:,ll),Nh,1)./(mpc.cp(1).*Uopt0(:,ll))  ).^(1/3);
    mpc.V      = reshape(temp,N,Nh); % rotor-averaged wind speed in the horizons

end

%% Assign the decision variables
Yopt          = Y;
Uopt          = Yopt(mpc.MU);
temp          = reshape(Uopt,[N,Nh]); % measured wind used
CT_prime      = temp(:,1);              % first action horizon
phi           = 0*[sol.turbine.Phi(1); zeros(N-1,1)];

