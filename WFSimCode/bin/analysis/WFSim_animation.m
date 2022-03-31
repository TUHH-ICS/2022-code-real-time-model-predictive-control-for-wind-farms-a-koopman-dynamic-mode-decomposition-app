function [hfig] = WFSim_animation( Wp,sol,hfig )
    % Import local variables from large structs
    Dr     = Wp.turbine.Drotor;
    ldyy   = Wp.mesh.ldyy;
    ldxx2  = Wp.mesh.ldxx2;
    yline  = Wp.mesh.yline;

    u_Inf  = Wp.site.u_Inf;

    N      = Wp.turbine.N;
    Cry    = Wp.turbine.Cry;
    Crx    = Wp.turbine.Crx;
    input  = sol.turbInput;

    time   = sol.time;
    k      = sol.k;
    
    turb_coord = .5*Dr*exp(1i*input.phi*pi/180);  % Yaw angles

    %% Create figure window, if not yet available
    if nargin <= 2
        scrsz = get(0,'ScreenSize');
        hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
    end;
    
    set(0,'CurrentFigure',hfig);
    
    idxX = 5:size(sol.u,1)-4;
    idxY = 5:size(sol.u,2)-4;

    %% Plot u velocity flow component
    subplot(2,2,[1 3]);
    contourf(ldyy(1,idxY ),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.2),'Linecolor','none');  
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.04]);  hold all; colorbar;
    axis equal; axis tight;
   
    % Plot the turbines in the field
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
        str = strcat('$T_{',num2str(kk),'}$');
        text(Cry(kk)+80,Crx(kk),str,'FontSize',16,'interpreter','latex')
    end
    % text(-600,ldxx2(end,end),['$t=~$ ', num2str(time,'%.1f'), '[s]'],'FontSize',20,'interpreter','latex');
    xlabel('$y$ [m]','FontSize',20,'interpreter','latex')
    ylabel('$x$ [m]','FontSize',20,'interpreter','latex');
    title(['$t=~$', num2str(time,'%.1f'),': Long. wind $u$ [m/s]'],'FontSize',20,'interpreter','latex');
    hold off;

    %% Plot the v velocity flow component
    subplot(2,2,2);
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',min(sol.v(idxX,idxY),u_Inf*1.2),'Linecolor','none');  colormap(hot);   hold all
    colorbar;
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
    end
    axis equal; axis tight
    xlabel('$y$ [m]','FontSize',20,'interpreter','latex')
    ylabel('$x$ [m]','FontSize',20,'interpreter','latex');
    title('Lat. wind $v$ [m/s]','FontSize',20,'interpreter','latex')
    hold off;

    %% Wake mean centreline first turbine
    D_ind    = yline{1};
    up(:,k)  = mean(sol.u(idxX,D_ind),2);

    set(0,'CurrentFigure',hfig);
    subplot(2,2,4)
    plot(ldxx2(idxX,1)',up(:,k),'k','Linewidth',2);
    xlabel('$x$ [m]','FontSize',20,'interpreter','latex');ylabel('$U_c$ [m/s]','FontSize',20,'interpreter','latex');grid;
    ylim([0 u_Inf+1]);xlim([ldxx2(idxX(1),1) ldxx2(idxX(end),1)]);
    title('Mean centerline wind $v$ [m/s]','FontSize',20,'interpreter','latex')
    for kk=1:N
        vline(Crx(kk),'r--')
        vline(Crx(kk)+3*Dr,'b:','3D')
    end
    hold off;
    drawnow;
end