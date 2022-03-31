function [CT_prime,phi,mpc] = qLPVMPCcontroller_dU_quadprog(sol,Wp,mpc,K,measured)

% initialize controller
mpc = MPCinit_ctrl(sol,Wp,mpc);
codedir = mfilename('fullpath');
maindir = fileparts(fileparts(fileparts(fileparts(codedir))));
dirFig = fullfile(maindir,'CtFig');    %define main directory
if ~isfolder(dirFig)
    mkdir(dirFig)
end
% solve mpc
if sol.k>=1
    
    % xinit = [Fk Pk CT_prime] in paper.
    xinit         = zeros(mpc.nx*Wp.turbine.N,1);
    xinit(mpc.Mf) = sol.turbine.force;
    xinit(mpc.Mp) = sol.turbine.power;
    xinit(mpc.Mu) = sol.turbine.CT_prime;
    
    % nl-1 is number of times the rotor-averaged wind speeds in the horizon
    % will be updated during one sample. If nl=1, the rotor-averaged wind speeds
    % are taken constant in the horizon
    nl = 1;
    Uopt = NaN(mpc.Nh*Wp.turbine.N,nl);
    
    for nlIdx = 1:nl
        
        % build wind farm model
        [K,mpc]  = wfmodel(sol,Wp,mpc,2,K);%wfmodel(sol,Wp,mpc,nlIdx,K);% nlIdx = 2 to estimate wind for Nh
        
        % build matrices horizon
        mpc  = Matsys(Wp,mpc,mpc.Nh);
        
        % solve the optimization problem and calculate time
        tic
        mpcest = mpc;
        mpcest.B = mpcest.Best;
        mpcest.BBt =mpcest.BBtest;
        mpcest = quadprog_sol(mpcest,Wp,sol,xinit);
        mpc = quadprog_sol(mpc,Wp,sol,xinit);
        mpc.displaytime = toc;
        
        %% True wind state space model of the wind farm
        aValue_approx = mpc.uval + mpc.utemp;
        X    = mpc.AA*xinit + mpc.BBt*aValue_approx;
        Y    = mpc.CC*X;                                %yopt
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh);  %Power output
        
        % tracking error value
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        
        % change in input error value
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpc.costValue = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;
        
        
        %% Estimated wind state space model of the wind farm with estims
        aValue_approx_est = mpcest.uval + mpcest.utemp;
        X_est    = mpcest.AA*xinit + mpcest.BBt*aValue_approx_est;
        Y_est    = mpcest.CC*X_est;                                %yopt
        P_est    = reshape(Y(mpcest.MP),Wp.turbine.N,mpcest.Nh);  %Power output
        
        % tracking error value
        Etemp_L    = mpcest.Pref(sol.k:sol.k+mpcest.Nh-1) - sum(P)';
        
        % change in input error value
        d_utemp2   = [aValue_approx(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approx(Wp.turbine.N+1:end) - aValue_approx(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpcest.costValue = Etemp_L'*mpcest.Q*Etemp_L + d_utemp2'*mpcest.R*d_utemp2;
        
        
    end
    
end
%% Assign the decision variables
Yopt          = Y;
Uopt          = Yopt(mpc.MU);
Yopt_est          = Y_est;
Uopt_est          = Yopt_est(mpcest.MU);

% if mod((sol.k-1)/100,1) == 0 %Figure for debugging
%     sN = 1;
%     figure;
%     for idxPl = 1:Wp.turbine.N, subplot(Wp.turbine.N,1,idxPl);
%         plot(Uopt(idxPl:Wp.turbine.N:end)), hold on;
%         plot(Uopt_est(idxPl:Wp.turbine.N:end),'--'); axis tight; grid on;
%         if idxPl == 1
%             title(sprintf('Time step %d: WT %d', sol.k,idxPl))
%         else
%             title(sprintf('WT %d', idxPl))
%         end
%         
%         if mod(idxPl-1,sN) == 0
%             ylabel('c_T [-]');
%         end
%         if idxPl > sN*(sN-1)
%             xlabel('k [-]');
%         end
%         if idxPl == Wp.turbine.N
%             legend('Opt','LCP');
%         end
%     end
%     print(gcf,fullfile(dirFig,sprintf('CTtrajectory%03d',sol.k)), '-dpng');
% end


if measured ==1
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]); % measured wind used
else
    temp          = reshape(Uopt_est,[Wp.turbine.N,mpc.Nh]);%estimated wind used
end
CT_prime      = temp(:,1);              % first action horizon
phi           = zeros(Wp.turbine.N,1);

%     figure;
%     for idxPl = 1:9,subplot(3,3,idxPl); plot(tempplot(idxPl:9:end)),...
%                 hold on;
%                 plot(aValue_approxLemke(idxPl:9:end),'--'); axis tight; grid on;...
%                 if idxPl == 1
%                     title(sprintf('Time step %d', sol.k))
%                 end
%     end

end
