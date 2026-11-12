function MPC_animation(sol_array,mpc)

%After running WFSim with the MPC, you can run this function to analyse
% the results

time = cell2mat(arrayfun(@(x)(x.k),sol_array,'UniformOutput', false));
te = 600;
Power = cell2mat(arrayfun(@(x)(x.turbine.power),sol_array,'UniformOutput', false));
N = size(Power,1);
CT_prime = cell2mat(arrayfun(@(x)(x.turbine.CT_prime),sol_array,'UniformOutput', false));

Pref_k = mpc.Pref(1:time(end));
sumPower = sum(Power); %e = Pref_k  - sumPower;

figure;
lw = 1; %LineWidth
cl = [lines(7);[0,0,0]; 0.5*ones(1,3)];
sl = repmat({'-';'--';'-.'},3,1);

subplot(2,1,1)
stairs(time,sumPower/1e6,'linewidth',lw);hold on;grid
stairs(time,Pref_k/1e6,'r--','linewidth',lw);
axis tight;
ylim([.8*min(Pref_k)/1e6 1.2*max(Pref_k)/1e6])
xlabel('$k$','interpreter','latex')
ylabel('$\sum P_i$ [MW]','interpreter','latex');
title('Wind farm power (blue) and reference (red)','interpreter','latex',...
'Fontsize',14);

subplot(2,1,2)
for kk = 1:N
    stairs(time,CT_prime(kk,:),'linewidth',lw,'Color',cl(kk,:),...
        'linestyle',sl{kk}); 
    hold on;
end
lc = regexp(sprintf('WT%d;',1:N),';','split');
legend(lc(1:end-1))
grid; axis tight; ylim([0 mpc.uM]) 
xlabel('$k$','interpreter','latex')
ylabel('$C_T^\prime$','interpreter','latex');

