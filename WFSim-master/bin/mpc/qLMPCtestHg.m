function [CT_prime,phi,mpc] = qLMPCtestHg(sol,Wp,mpc)
% qLMPCcontroller calculates optimal MPC CT_prime values based on current
% effective wind speeds and thrust coeffecients via QP

dirFig = fullfile(sol.OutDir,'CtFigH');    %define main directory
if ~isfolder(dirFig)
    mkdir(dirFig)
end

%% Parameter for optimization
Nh = mpc.Nh;        % prediction trajectory length
N = Wp.turbine.N;   % number turbines

% QP boundary condition: (Aco_) x <= bco_
Ac0_ = [eye(Nh*N); -eye(Nh*N)];
bc0_ = [mpc.Ulim_lower; - mpc.Ulim_upper];

% xinit = [Fk Pk CT_prime] 
xinit         = zeros(mpc.nx*Wp.turbine.N,1);
xinit(mpc.Mf) = sol.turbine.force;
xinit(mpc.Mp) = sol.turbine.power;
xinit(mpc.Mu) = sol.turbine.CT_prime;

% nl-1 is number of times the rotor-averaged wind speeds in the horizon
% will be updated during one sample. If nl = 1, the rotor-averaged wind speeds
% are taken constant in the horizon
nl = 1;
Uopt0 = NaN(Nh*N,nl);
Uopt1 = NaN(Nh*N,nl);
%DT = Wp.sim.h;
mpc1 = mpc;

for ll = 1:nl
    
    %% Default mpc as QP:  Build wind farm model
    mpc  = wfmodel(sol,Wp,mpc,ll);
    
    % build matrices horizon, Hessian, and gradient 
    mpc  = matsysLSHg(Wp,mpc,sol,xinit);
    % solve QP optimization problem and calculate time
    tic
    bc   = - bc0_ + Ac0_ * mpc.utemp;
    [uval,fval] = myquadprog(mpc.H,mpc.g,-Ac0_,bc);
    mpc.x = uval + mpc.utemp;
    mpc.uval = uval;
    mpc.fval = fval;
    mpc.displaytime = toc;
    
    % True wind state space model of the wind farm
    uOpt = mpc.uval + mpc.utemp;
    X    = mpc.AA*xinit + mpc.BBt*uOpt;
    Y    = mpc.CC*X;                                %yopt
    P    = reshape(Y(mpc.MP),N,Nh);  %Power output
    
    % tracking error value
    Etemp_L    = mpc.Pref(sol.k:sol.k+Nh-1) - sum(P)';
    
    % change in input error value
    d_utemp2 = [uOpt(1:N) - sol.turbine.CT_prime;...
        uOpt(N+1:end) - uOpt(1:end-N)];
    
    % Calculate the cost function
    mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
    
    Uopt0(:,ll) = uOpt; %value(Y(mpc.MU))
    Popt       = P;
    
    temp       = (repmat(Popt(:,ll),Nh,1)./(mpc.cp(1).*Uopt0(:,ll))  ).^(1/3);
    mpc.V      = reshape(temp,N,Nh); % rotor-averaged wind speed in the horizons
    
    %% HSqLd_WF: build wind farm matrices A and B, stacked prediction horizon L and S, 
    % and Hessian H, gradient g
%   mpc2 = HSqLd_WF(sol,Wp,mpc2,ll,xinit);
    mpc1 = HSqLdD_WF(sol,Wp,mpc1,ll,xinit,mpc);
    % mpc1 = HSqLdPhi_WF(sol,Wp,mpc1,ll,xinit);
    
    % solve QP optimization problem and calculate time
    tic
    bc   = - bc0_ + Ac0_ * mpc1.utemp;
    %[uval2,fval2] = myquadprog(mpc1.H,mpc1.g,-Ac0_,bc);
    [uTilde,fval1] = myquadprog(mpc1.H,mpc1.g,-Ac0_,bc);
    mpc1.x = uTilde + mpc1.utemp;
    mpc1.uval = uTilde;
    mpc1.fval = fval1;
    mpc1.displaytime = toc;
    
    uOpt = uTilde + mpc.utemp; % utilde = uOpt - utemp: Get ctrl input of orig problem
    
    % True wind state space model of the wind farm
    nxWF = mpc.nx *N;
    mpc_MPV = repmat([mpc.MP(1:nxWF); zeros(nxWF,1)],Nh,1);
    mpc_MUV = repmat([mpc.MU(1:nxWF); zeros(nxWF,1)],Nh,1);
       
    % change in input error value
    d_utemp2 = uOpt - [sol.turbine.CT_prime; uOpt(1:(N*(Nh-1)))]; %delta U from U
    X = mpc1.AA*mpc1.xinitD + mpc1.BBt*d_utemp2;
    Y1 = X;                                %yopt
    P1 = reshape(Y1(mpc_MPV == 1),N,Nh);  %Power output
    
    % tracking error value
    Etemp_L    = mpc1.Pref(sol.k:sol.k+Nh-1) - sum(P1)';
    
    % Calculate the cost function
    mpc1.costValue = Etemp_L'*mpc1.Q*Etemp_L + d_utemp2'*mpc1.R*d_utemp2;
    
    Uopt1(:,ll) = uOpt; %Y1(mpc_MUV==1);
    Popt       = P1;
    
    temp       = (repmat(Popt(:,ll),Nh,1)./(mpc1.cp(1).*Uopt1(:,ll))  ).^(1/3);
    mpc1.V      = reshape(temp,N,Nh); % rotor-averaged wind speed in the horizons

end

%% Assign the decision variables
Yopt          = Y;
Uopt          = Yopt(mpc1.MU);
temp          = reshape(Uopt,[N,Nh]); % measured wind used
CT_prime      = temp(:,1);            % first action horizon
phi           = zeros(N,1);

if mod((sol.k-1)/100,1) == 0 %Figure for debugging
    sN = 1;
    figure;
    for idxPl = 1:N, subplot(3,3,idxPl); %N
        plot(Uopt0(idxPl:N:end,nl)), hold on;
        plot(Uopt1(idxPl:N:end,nl),'--'); axis tight; grid on;
        if idxPl == 1
            title(sprintf('Time step %d: WT %d', sol.k,idxPl))
        else
            title(sprintf('WT %d', idxPl))
        end

        if mod(idxPl-1,sN) == 0
            ylabel('c_T [-]');
        end
        if idxPl > sN*(sN-1)
            xlabel('k [-]');
        end
        if idxPl == N
            legend('QP','velQP');
        end
    end
    print(gcf,fullfile(dirFig,sprintf('CTtrajectory%03d',sol.k)), '-dpng');
end


