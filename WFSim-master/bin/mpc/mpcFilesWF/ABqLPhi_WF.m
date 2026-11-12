function [Aeval,Beval] = ABqLPhi_WF(rhoP,Phi,mpc)
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
% modified for WF simulation: Antje Dittmer


%% Parameters %if length(rhoP) == 1
N = size(Phi,1);
nx = 3; 
Vbar = rhoP(:,1); % Vbar

%% Tr=0.5*pi*rb^3*Cq*v^2
% controller models
tau  = mpc.tau; % 5;                             % time constant filter CT'

%% Create system matrix A:
% WT: States: [Theta_s [rad], omega_r [rad/s], omega_gr [rad/s], T_g [Nm]]
% WT: States for model xA = [y, xdot] with y = [Tg, beta]
% States: [F, P, CT'];
% WT: States for model xA = [y, xdot] with y = [x] = [F, P, CT']
Atemp = -1/tau * eye(N * nx);

%% Create input matrix B
% Input orig model: V [m/s], beta [rad], T_g_ref [Nm]
% -> vel model: Vdot [m/s^2]  betaDot [rad/s], T_g_refDot [Nm/s]
% for beta: Vdot [m/s^2]  T_g_refDot [Nm/s] betaDot [rad/s], 
%Btemp = [zeros(1,3); dTrdV/Jr, 0 dTrdbeta/Jr; zeros(1,3);  0,  1/tau_g,0];
Bvec = 1/tau * [- mpc.cf .* (cos(Phi).* Vbar).^2, mpc.cp .* (cos(Phi).*Vbar).^3 ,ones(N,1)]'; %changing init with zeros
Btemp = zeros(N*nx,N);
for idx = 1:N
    Btemp((idx-1)*nx+(1:nx),idx)= Bvec(:,idx);
end
  
% With Beta
Aeval = Atemp;
Beval = Btemp; 

end

%% Physical and model constants
% rho = 1; %1.2250;
% rb = 63; % m blade radius
% Ng = 97; %Gearbox ratio
% Jr = (115926 + 3 * 11.776e6) *1; %kg*m^2 Inertia of the rotor (Hub inertia + 3 blades
% Jg = 534.116;% * Ng^2; %kg*m^2 Inertia of the generator (HSS)
% Jgr = Jg * Ng^2; % Inertia of the generator in LSS
% Ks = 867637000; %Nm/rad Stiffness of the transmission
% Bs = 6215000;  %Nm/rad/sec        %Damping of the transmission
% tau_g = 0.01; % for pitch actuator and airspeed
