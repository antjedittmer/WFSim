function [Aeval,Beval] = ABqL_WT_beta(rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ)
% ABqL calculates system (A) and input (B) matrix for quasiLPV
% representation for ve
%
% Inputs: * Ts: Sampling time
%         * AwithV: Scheduling parameter vector: AwithV = V
%
% Outputs: * A: System matrix A(AwithV) at current estimated air speed AwithV, discretised with Ts
%          * B: Input matrix B(AwithV) at current AwithV, discretised with Ts

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WECS simulation: Antje Dittmer

persistent lambdaVec pitchVec cQ Cqdlambda Cqdbeta

if isempty(lambdaVec)
    lambdaVec =  Rotor_Lambda;
    pitchVec = Rotor_Pitch;
    cQ = Rotor_cQ;
    %% Compute gradients %%[Cpdlambda,Cpdbeta] = gradient(Cp,mean(diff(Xq(1,:))),mean(diff(Yq(:,1))));
    [Cqdlambda,Cqdbeta] = gradient(Rotor_cQ,mean(diff(Rotor_Lambda)),mean(diff(Rotor_Pitch))); %dcQ/dbetaDeg * dbetaDeg/dbetaRad *betaRAd
end
   

%% Parameters
if length(rhoP) == 1 
    Vbar = rhoP; % Vbar
    LUT_V = 4:25; % LUT table rows
    %LUT_omega = [0.4815, 0.6017,0.7217,0.8413,0.9603,1.0786,1.1949, repmat(1.2671,1,15)];
    LUT_omega =  [0.7459    0.7776    0.8207    0.8738    0.9354    1.0429    1.1592    1.2371 ....
         1.2554    1.2548    1.2545    1.2546    1.2545    1.2545    1.2545    1.2545    1.2545    1.2545 ...
    1.2545    1.2545    1.2545    1.2545];
    
%     LUT_beta = [zeros(1,8),0.0642    0.1126    0.1493    0.1808    0.2093    0.2356    0.2602...
%         0.2833    0.3054    0.3269    0.3477    0.3676    0.3868    0.4056];  
%     
    LUT_beta = [zeros(1,8), 0.0642    0.1126    0.1493    0.1808    0.2093    0.2370    0.2620    0.2854    0.3077    0.3293 ...
    0.3503    0.3704    0.3897    0.4086];
    
    % Get omegabar and beta
    omegabar = interp1(LUT_V,LUT_omega,Vbar,'pchip');
    beta_bar = interp1(LUT_V,LUT_beta,Vbar,'pchip');
    
elseif length(rhoP) == 3
    % rho = xdotk_1([1,3,4]); % Tg not included: [Omega [1/s], beta[rad], V [m/s]];
    Vbar = rhoP(3);
    omegabar = rhoP(1); % pass OmegaR OmegaG/Ng as approximation of omegaR
    beta_bar = rhoP(2);
    
else
    error('Unexpected length of scheduling parameter rho')
end

% LUT_torque = [2.2068    4.9984    8.8210   13.5123   18.9661   23.8577,...
%     29.4753   37.0202   43.5290   43.5290   43.5290   43.5290   43.5290,...
%     43.5300   43.5300   43.5300   43.5300   43.5300,...
%    43.5300   43.5300   43.5300   43.530];

%% Physical and model constants
rho = 1.2250;
rb = 63; % m blade radius
Ng = 97; %Gearbox ratio
Jr = (115926 + 3 * 11.776e6) *1; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
Jgr = Jg * Ng^2; % Inertia of the generator in LSS
Ks = 867637000; %Nm/rad Stiffness of the transmission
Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
tau_g = 0.01; % for pitch actuator and airspeed

%% Lambda 
Lambda_bar = omegabar*rb/Vbar;
dlambdadomega= rb/Vbar;
dlambdadV = - omegabar*rb/Vbar^2;

%% Grid points
Cqb = interp2(lambdaVec,pitchVec,cQ,Lambda_bar,beta_bar); 
dCqdlambdab = interp2(lambdaVec,pitchVec,Cqdlambda,Lambda_bar,beta_bar);
dCqdbetab = interp2(lambdaVec,pitchVec,Cqdbeta,Lambda_bar,beta_bar);

%% Tr=0.5*pi*rb^3*Cq*v^2
dTrdomega = 0.5*rho*pi*rb^3*(dCqdlambdab*dlambdadomega)*Vbar^2;% br = dTr/domega
dTrdbeta = 0.5*rho*pi*rb^3*dCqdbetab*Vbar^2;% Krb = dTr/dbeta
dTrdV =  2*0.5*rho*pi*rb^3*Vbar*Cqb + ...
    0.5*rho*pi*rb^3*Vbar^2*dCqdlambdab*dlambdadV;

%% Create system matrix A:
% States: [Theta_s [rad], omega_r [rad/s], omega_gr [rad/s], T_g [Nm]]
% States for model xA = [y, xdot] with y = [Tg, beta]
% x = [ThetaS, OmegaR, OmegaG, Tg, beta,]
Atemp = [[0,1,-1,0] ;
    [-Ks/Jr, (-Bs + dTrdomega)/Jr, Bs/Jr,0];
    Ks/Jgr, Bs/Jgr, -Bs/Jgr, - Ng/Jgr; 
    zeros(1,3), -1/tau_g];

%% Create input matrix B
% Input orig model: V [m/s], beta [rad], T_g_ref [Nm]
% -> vel model: Vdot [m/s^2]  betaDot [rad/s], T_g_refDot [Nm/s]
% for beta: Vdot [m/s^2]  T_g_refDot [Nm/s] betaDot [rad/s], 
Btemp = [zeros(1,3); dTrdV/Jr, 0 dTrdbeta/Jr; zeros(1,3);  0,  1/tau_g,0];

% torque_bar = interp1(4:25,LUT_torque,Vbar,'pchip');
% u0 = [Vbar; torque_bar*10^3; beta_bar;];
% x0 = [0,omegabar, omegabar,  u0(2)]';
% 
% tmp = (-Btemp(2,:)*u0 - Atemp(2,[1,3])* x0([1,3]))/omegabar,
% Atemp(2,2) = tmp;
% x0_sys5DoF = - inv(Atemp) * Btemp * u0;


% With Beta
Aeval = [Atemp,Btemp(:,3); zeros(1, size(Atemp,2)), -1/tau_g];
Beval = [Btemp; zeros(1, size(Btemp,2))];
Beval(:,3) = 0;
Beval(end,3) = 1/tau_g;

end
