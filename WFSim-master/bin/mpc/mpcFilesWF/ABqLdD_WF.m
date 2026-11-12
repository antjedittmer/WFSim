function [Av,Bv] = ABqLdD_WF(DT,rhoP,Phi,mpc)
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

% Calculate system and input matrices at current wind speed
Vbar = rhoP(:,1);
[Arho,Brho] = ABqL_WF(Vbar,0*Phi,mpc);

%% Calculate discrete forms of A,B
A2 = Arho*Arho;
A3 = Arho*A2;
Apre = (1.0/24.0)*DT^4*A3+(1.0/6.0)*DT^3*A2+0.5*DT^2*Arho+DT*eye(size(Arho));

%% Finish calculation of A,B
Ak = Apre*Arho + eye(size(Arho));
Bk = Apre*Brho;
% sysd = c2d(ss(Arho,Brho,eye(size(Arho)), zeros(size(Arho,1), size(Brho,2))),DT);

%% Calculate new matrices augmented system
% Av = Ak; %sysd.A; % 
% Bv = Bk; % sysd.B; %

ny = mpc.N*3; %
nx = size(Arho,1);
% nyvec = 2:3: mpc.N*3; % 
% nd = 1;% B: columns disturbance input
% nuvec = 2:3; % B: columns control input
Cdiag = eye(nx);%kron(eye(mpc.N),[0 1 0]);

Av = [eye(ny), Cdiag * Ak;
    zeros(size(Ak,1),ny), Ak];
Bv = [Cdiag * Bk; Bk];
% Bv1 = [Cdiag * Bk(nyvec,nd); Bk(:,nd)];

% Av = [A,Bv1; [zeros(nd,(nx+ny)) 0*eye(nd)]];
% Bv = [B; zeros(nd, size(B,2))];


end