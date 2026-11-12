function [Av,Bv] = ABqLd_WF(DT,rhoP,Phi,mpc,isVelForm)
% ABqLd calculates system (A) and input (B) matrix for quasiLPV
% representation from ABgL and augments the states with integrating
% disturbances.
%
% Inputs: * DT: Sampling time
%         * rho: Scheduling parameter vector: rho = V
%
% Outputs: * A: Augmented system matrix A(rho) at current rho
%          * B: Augmented Input matrix B(rho) at current rho

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WF simulation: Antje Dittmer

%% Set standard or velocity form
if nargin < 5
    isVelForm = 0;
end

%% Calculate continous system and input matrices at current wind speed
Vbar = rhoP(:,1);
[Ac,Bc] = ABqL_WF(Vbar,0*Phi,mpc);

if  isVelForm == 0
    Arho = Ac;
    Brho = Bc;
    
else
    ny = size(Ac,1); %nx = ny
    nu = size(Bc,2); 
    Arho = [zeros(ny),eye(ny); zeros(ny), Ac];
    Brho = [zeros(ny,nu); Bc];
    
end

%% Calculate discrete forms of A,B
A2 = Arho*Arho;
A3 = Arho*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*Arho+DT*eye(size(Arho));

%% Finish calculation of A,B
Av = Apre*Arho + eye(size(Arho));
Bv = Apre*Brho;
% sysd = c2d(ss(Arho,Brho,eye(size(Arho)), zeros(size(Arho,1), size(Brho,2))),DT);
% Av1 = sysd.A;
% Bv1 = sysd.B;
%[Av,Bv,Av1,Bv1] = ABqLd_WF(DT,rhoP,Phi,mpc,isVelForm)
end