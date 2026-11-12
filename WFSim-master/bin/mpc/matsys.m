function mpc = matsys(Wp,mpc)
% matsys sets up the stacked matrices AA/L, BBt/S and CC

%% Assign inputs
A = mpc.A;
Bt = mpc.Bt;
C  = mpc.C;

Nh = mpc.Nh;
nx = size(A,1);
N = Wp.turbine.N;

%% Calculate stacked matrices L and S
L = zeros(nx*Nh,nx);
S = zeros(nx*Nh,N*Nh);

L(1:nx,1:nx) = A;
S(1:nx,1:N) = Bt(:,:,1);
for idx = 1:Nh-1
    L(idx*nx+1:(idx+1)*nx,1:nx) = A*L((idx-1)*nx+1:idx*nx,1:nx);
    S(idx*nx+1:(idx+1)*nx,1:(idx+1)*N) = [A*S((idx-1)*nx+1:idx*nx,1:idx*N),Bt(:,:,idx+1)];
end

Cc = kron(eye(Nh),C);

%% Assign outputs
mpc.AA  = L;
mpc.BBt = S;
mpc.CC  = Cc;

end

