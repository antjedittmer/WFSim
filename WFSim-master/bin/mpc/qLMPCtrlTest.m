function [CT_prime,phi,mpc] = qLMPCtrlTest(sol,Wp,mpc)

% Controller and dirFig for comparison with orig. controller
mpc1 = mpc; % initialize controller for NMPC

dirFig = fullfile(sol.OutDir,'CtFig');    %define main directory
if ~isfolder(dirFig)
    mkdir(dirFig)
end

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
Uopt1 = NaN(Nh*N,nl);

for ll = 1:nl
    
    % build wind farm model
    mpc  = wfmodel(sol,Wp,mpc,ll);
    
    % build matrices horizon, Hessian, and gradient 
    mpc  = matsysLSHg(Wp,mpc,sol,xinit);
    %sol.turbInput.phi = 0 * sol.turbInput.phi;
    %[H,g,mpc] = HSqLd_WF(sol,Wp,mpc,ll,xinit);
    
    % solve QP optimization problem and calculate time A(u_tilde +
    tic
    bc   = bc0_ - Ac0_ * mpc.utemp; % A(u_tilde + utemp) <= Bco
    [uval,fval] = myquadprog(mpc.H,mpc.g,Ac0_,bc);
    [uval0,fval0] = myquadprog(mpc.H,mpc.g0,Ac0_,bc0_);

    mpc.QPerr = max(abs(uval0 - (uval + mpc.utemp)));
    mpc.x = uval0; % uval + mpc.utemp;
    mpc.uval = uval;
    mpc.fval = fval;
    mpc.fval = fval0;
    mpc.displaytime = toc;
    
    %% True wind state space model of the wind farm
    uTilde = mpc.uval + mpc.utemp;
    X    = mpc.AA*xinit + mpc.BBt*uTilde;
    Y    = mpc.CC*X;                                %yopt
    P    = reshape(Y(mpc.MP),N,Nh);  %Power output
    
    % tracking error value
    Etemp_L = mpc.Pref(sol.k:sol.k+Nh-1) - sum(P)';
    
    % change in input error value
    d_utemp2 = [uTilde(1:N) - sol.turbine.CT_prime;...
        uTilde(N+1:end) - uTilde(1:end-N)];
    
    % Calculate the cost function
    mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
    
    Uopt0(:,ll) = uTilde; %value(Y(mpc.MU))
    Popt       = P;
    
    temp       = (repmat(Popt(:,ll),Nh,1)./(mpc.cp(1).*Uopt0(:,ll))  ).^(1/3);
    mpc.V      = reshape(temp,N,Nh); % rotor-averaged wind speed in the horizons

    %% Optimal trajectory calculated with optimize
    yalmip('clear');
    
    % define decision variables for the windfarm
    U    = sdpvar(N*Nh,1);
    
    % build wind farm model
    mpc1  = wfmodel(sol,Wp,mpc1,ll);
    
    % build matrices horizon
    mpc1  = matsys(Wp,mpc1);
    
    X    = mpc1.AA*xinit + mpc1.BBt*U ;
    Y1   = mpc.CC*X;
    P1   = reshape(Y1(mpc.MP),N,Nh) ;
    
    E    = mpc1.Pref(sol.k:sol.k+Nh-1)-sum(P1)';
    
    cons = mpc.um <= U <= mpc.uM;
    
    dU   = [ U(1:N)-sol.turbine.CT_prime ; U(N+1:end)-U(1:end-N)];
    
    cost = E'*mpc.Q*E + dU'*mpc.R*dU;
    
    ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1);
    
    optimize(cons,cost,ops);
    
    Uopt1(:,ll) = value(U); %value(Y(mpc.MU))
    Popt       = value(P1);
    
    temp       = ( repmat(Popt(:,ll),Nh,1)./(mpc.cp(1).*Uopt1(:,ll))  ).^(1/3);
    mpc1.V      = reshape(temp,N,Nh); % rotor-averaged wind speed in the horizons

end

%% Assign the decision variables
Yopt          = Y;
Uopt          = Yopt(mpc.MU);
temp          = reshape(Uopt,[N,Nh]); % measured wind used
CT_prime      = temp(:,1);              % first action horizon
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
            legend('QP','Opt');
        end
    end
    print(gcf,fullfile(dirFig,sprintf('CTtrajectory%03d',sol.k)), '-dpng');
end



%     figure;
%     for idxPl = 1:9,subplot(3,3,idxPl); plot(tempplot(idxPl:9:end)),...
%                 hold on;
%                 plot(aValue_approxLemke(idxPl:9:end),'--'); axis tight; grid on;...
%                 if idxPl == 1
%                     title(sprintf('Time step %d', sol.k))
%                 end
%     end

end
