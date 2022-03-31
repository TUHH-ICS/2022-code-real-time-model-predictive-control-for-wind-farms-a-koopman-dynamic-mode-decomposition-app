function [JR,JQ] = plotWFSimTimeseries(Wp,sol_array,controller,Power,CT_prime,Phi,V,JR,time,measured,figName,dirFig,mpc,K)

if nargin < 13
    controller = 0;
    mpc ='';
end

if nargin < 14
    K = [];
end

JQ = nan;
N = Wp.turbine.N;
fs = 12; % Frontsize for plots;
% Assign variables from sol_array structure
kk = length(Power);

% Values closed loop
if controller > 0
    P_greedy = mpc.Pgreedy*ones(kk,1); % 2.876e+06in MW
    e = mpc.Pref(1:kk) - sum(Power)';
    tVAF = 15;
    JQ = mean(abs(e(tVAF:end)));
end

if controller > 1
    X_est = cell2mat(arrayfun(@(x)(x.xprev),sol_array,'UniformOutput', false)); %X_est = cell2mat(arrayfun(@(x)(x.x),sol_array,'UniformOutput', false));
    V_est = X_est(1:2,:);
    VAF1 = 1-var(V(1,tVAF:end)-V_est(1,tVAF:end))/var(V(1,tVAF:end));
    VAF2 = 1-var(V(2,tVAF:end)-V_est(2,tVAF:end))/var(V(2,tVAF:end));
end

figure; 
% Prepare lines etc.
% posDefault = get(2, 'position');
posDefault = [805.0000  272.2000  560.0000  420.0000];
set(gcf, 'position', [posDefault(1:2),posDefault(3)*1.2,posDefault(4)*1.2]);
cl = colormap('lines');
nl = length(unique(cl(:,1)));
cl(nl+1,:) = 0.5 *ones(1,3);
cl(nl+2,:) = 0 *ones(1,3);
%fsl = 10;
fsL = 10;
lw = 1; % linewidth
ls = repmat({'-','--','-.'},1,3);
if controller == 0
    titlestring = sprintf('OL %s AA=%2.2e',figName, JR); 
elseif controller == 1 || measured == 1
    %     filename = sprintf('V_inf=%.0f: Q=%0.1d R=%g T.E=%2.2e[W] A.A=%2.2e',...
    %         Wp.site.u_Inf,mpc.Q(1,1),mpc.R(1,1),JQ,JR);
    titlestring = sprintf('Q=%0.1d, R=%g: TE=%2.2e[W] AA=%2.2e',...
        mpc.Q(1,1),mpc.R(1,1),JQ,JR);
else
    %     filename = sprintf('V_inf=%.0f States=%d: Q=%0.1d R=%g T.E=%2.2e[W] A.A=%2.2e',...
    %         Wp.site.u_Inf,length(K.A),mpc.Q(1,1),mpc.R(1,1),JQ,JR);
    titlestring = sprintf('$no_{States}$=%d, Q=%0.1d, R=%g: TE=%2.2e[W] AA=%2.2e',...
        length(K.A),mpc.Q(1,1),mpc.R(1,1),JQ,JR);
end

subplot(4,1,1)
if controller == 0
    stairs(time,sum(Power)/1e6,'--','linewidth',lw+0.2,'color',cl(4,:,:));
    axis tight; grid on;
    legend('$P_T$', 'Orientation','horizontal','interpreter','latex',...
        'Fontsize',fsL)
else
    stairs(time,mpc.Pref(1:kk)/1e6,'k','linewidth',lw);
    hold on;grid on;
    stairs(time,sum(Power)/1e6,'--','linewidth',lw+0.2,'color',cl(4,:,:));
    hold on;grid on;
    stairs(time,P_greedy(1:kk)/1e6,'--','linewidth',lw+0.2,'color',cl(5,:,:));
    axis tight;
    ylim([.8*min(mpc.Pref(1:kk))/1e6 1.2*max(mpc.Pref(1:kk))/1e6])
    legend('$P_{ref}$','$P_T$','$P_{greeedy}$', 'Orientation','horizontal','interpreter','latex',...
        'Fontsize',fsL)
end
ylabel('$\sum P_i$ [MW]','interpreter','latex','Fontsize',fs);
title(titlestring, 'Fontsize',fs,'interpreter','latex')

subplot(4,1,2)
for idx = 1:N
    stairs(time,CT_prime(idx,:),'-','color',cl(idx,:),'linewidth',lw);hold on;
end
grid on; axis tight; ylim([0.1 2.1])% 
ylabel('Input $C_T^\prime$ [-]','interpreter','latex','Fontsize',fs);
%title({'\bf First Turbine (blue) and Second Turbine (red)',['Actuator activity= ',num2str(JR)]});
strLeg = sprintf('WT%d,',1:N);
legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fsL);

subplot(4,1,3)
usePhi = mean(diff(Phi(1,6:end))) > 0;
if usePhi
    for idx = 1:N
        stairs(time,Phi(idx,:),ls{idx},'color',cl(idx,:),'linewidth',lw);hold on;
        ylabel('Input $\Phi$ [deg]','interpreter','latex','Fontsize',fs);
        strLeg = sprintf('WT%d,',1:N);
        legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fsL);
    end
elseif measured == 1 || controller == 0 || controller == 1
    stairs(time,V(1,:),ls{1},'color',cl(1,:),'linewidth',lw);hold on; grid on
    ylabel('Wind $U_{r}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT1','interpreter','latex','Fontsize',fsL);
else
    stairs(time,V(1,:),ls{1},'color',cl(1,:),'linewidth',lw);hold on;
    stairs(time,V_est(1,:),'--','linewidth',lw, 'color',[0,0,0.7]);hold on;
    
    ylabel('Wind $U_{r}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT1',['WT1 est.,VAF:',num2str(round(VAF1*100)),'$\%$'],...
        'Orientation','horizontal',...
        'interpreter','latex','Fontsize',fsL,'Location','northeast');
    pos = axis;
    axis([pos(1:3),1.2*pos(4)])
    
end
grid on; axis tight;

subplot(4,1,4)
if (measured == 1 || controller == 0 || controller == 1) && usePhi
    stairs(time,V(1,:),ls{1},'color',cl(1,:),'linewidth',lw);hold on; grid on
    stairs(time,V(2,:),ls{1},'color',cl(2,:),'linewidth',lw);hold off; grid on
    xlabel('$k$ [-]','interpreter','latex','Fontsize',fs)
    ylabel('Wind $U_{r}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT1','WT2','interpreter','latex','Fontsize',fsL);
    
elseif usePhi
    stairs(time,V(1,:),ls{1},'color',cl(1,:),'linewidth',lw);hold on;
    stairs(time,V_est(1,:),'--','linewidth',lw, 'color',[0,0,0.7]);hold on;
    stairs(time,V(2,:),ls{1},'color',cl(2,:),'linewidth',lw);hold on;grid on
    stairs(time,V_est(2,:),'--','linewidth',lw,'color',[0.7,0,0]);hold off;
    
    xlabel('$k$ [-]','interpreter','latex','Fontsize',fs)
    ylabel('Wind $U_{r1}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT1',['WT1 est.,VAF:',num2str(round(VAF1*100)),'$\%$'],...
        'WT2',['WT2 est., VAF:',num2str(round(VAF2*100)),'$\%$'],...
        'Orientation','horizontal',...
        'interpreter','latex','Fontsize',fsL,'Location','northeast');
    
elseif (measured == 1 || controller == 0 || controller == 1)
    
    stairs(time,V(2,:),ls{1},'color',cl(2,:),'linewidth',lw);grid on
    
    xlabel('$k$ [-]','interpreter','latex','Fontsize',fs)
    ylabel('Wind $U_{r}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT2',...
        'Orientation','horizontal',...
        'interpreter','latex','Fontsize',fsL,'Location','northeast');
    
else
    stairs(time,V(2,:),ls{1},'color',cl(2,:),'linewidth',lw);hold on;grid on
    stairs(time,V_est(2,:),'--','linewidth',lw,'color',[0.7,0,0]);hold off;
    axis tight;
    
    xlabel('$k$ [-]','interpreter','latex','Fontsize',fs)
    ylabel('Wind $U_{r}$ [m/s]','interpreter','latex','Fontsize',fs);
    legend('WT2',['WT2 est., VAF:',num2str(round(VAF2*100)),'$\%$'],...
        'Orientation','horizontal',...
        'interpreter','latex','Fontsize',fsL,'Location','northeast');
end
axis tight;
pos = axis;

axis([pos(1:3),1.2*pos(4)])

%filenamepng = strrep(fullfile(dirFig,figName),'.','dot');
filenamepng = matlab.lang.makeValidName(figName);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');