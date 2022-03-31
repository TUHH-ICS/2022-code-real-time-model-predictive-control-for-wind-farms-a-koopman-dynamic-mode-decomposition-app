function fig1 = plotEDMDinputsEffWind(ysim_val,Inputs_val,Deterministic_val, FITje_val,dirFig)

    fig1 = figure;
    cl = colormap('lines');
    lw = 1;
    fs = 10;
    kplot = 60: length(ysim_val);
    
    noSub = size(Inputs_val,1) + 1;
    subplot(noSub,1,1)
    plot(Inputs_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    hold on;
    plot(Inputs_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)) %3
    grid on; axis tight; ylim([0 2+0.2])
    ylabel('Input $C_T$ [-]','interpreter','latex','Fontsize',fs);
    title(sprintf('Validation Set: Inputs and Effective Wind: VAF Ur1: %d, Ur2: %d', round(FITje_val')),'interpreter','latex','Fontsize',fs);
    strLeg = sprintf('WT%d,',1:2);
    legend(regexp(strLeg(1:end-1), ',', 'split'),'interpreter','latex','Fontsize',fs, 'Orientation','horizontal', 'location','southeast'); 
    legend('boxoff')
    set(gca,'fontsize', fs)
    
    if noSub > 3
            subplot(noSub,1,2)
            plot(Inputs_val(3,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
            grid on; axis tight
            ylabel('Input $\Phi$ [deg]','interpreter','latex','Fontsize',fs);
            set(gca,'fontsize', fs)
    end

    subplot(noSub,1,noSub-1)
    plot(Deterministic_val(1,kplot),'LineWidth',lw,'Color', cl(1,:)); %1
    hold on;
    plot(ysim_val(kplot,1),'Color',[0,0,0.6],'LineWidth',lw,'LineStyle','--') %3
    grid on; axis tight;
    pos = axis; axis([pos(1:3), pos(4)+0.3]);
    %xlabel('k [-]','interpreter','latex')
    ylabel('$U_{r1}$ [m/s]','interpreter','latex')
    legend({'Real','Est.'},'Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', fs)

    subplot(noSub,1,noSub)
    plot(Deterministic_val(2,kplot),'LineWidth',lw,'Color', cl(2,:)); %1
    hold on;
    plot(ysim_val(kplot,2),'Color',[0.6, 0, 0.0],'LineWidth',lw,'LineStyle','--') %3
    grid on; axis tight
    pos = axis; axis([pos(1:3), pos(4)+0.3]);
    xlabel('k [-]','interpreter','latex')
    ylabel('$U_{r2}$ [m/s]','interpreter','latex')
       
    legend({'Real','Est.'},'Location','northeast','interpreter','latex','fontsize',fs,'orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', fs)
    
    %sys_red{1,1}.Xint = xo;
    
    print(gcf,fullfile(dirFig,'Val'), '-dpng');
    print(gcf,fullfile(dirFig,'Val'), '-depsc');