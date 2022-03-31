function [CT_prime,phi,mpc] = qLPVMPCcontroller_test(sol,Wp,mpc)

% initialize controller
mpc = MPCinit_ctrl(sol,Wp,mpc);

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
        
        yalmip('clear');
              
        % define decision variables for the windfarm
        U    = sdpvar(Wp.turbine.N*mpc.Nh,1);
        
        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,nlIdx);
        
        % build matrices horizon
        mpc  = matsys(Wp,mpc);      
        
        % prepare in state space of wind farm
        X    = mpc.AA*xinit + mpc.BBt*U ;
        Y    = mpc.CC*X;
        
        %Power output
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh) ;
        
        % tracking error
        E    = mpc.Pref(sol.k:sol.k+mpc.Nh-1)-sum(P)';
        
        % cons was initialized as empty above
        cons = mpc.um <= U <= mpc.uM; 
        
        % change in input
        dU   = [ U(1:Wp.turbine.N)-sol.turbine.CT_prime ; U(Wp.turbine.N+1:end)-U(1:end-Wp.turbine.N)];%CT
        
        % cost function 
        cost = E'*mpc.Q*E + dU'*mpc.R*dU;
        
        % YALMIP toolbox using the sdpsetting and optimizer function
        tic;
        ops  = sdpsettings('solver','','verbose',0,'cachesolvers',1);%sdpsettings('solver','cplex','verbose',0,'cachesolvers',1)
        optimize(cons,cost,ops);
        tsolver = toc;
        Uopt(:,nlIdx) = value(U); %value(Y(mpc.MU))
        Popt       = value(P);
        temp       = (repmat(Popt(:,nlIdx),mpc.Nh,1)./(mpc.cp(1).*Uopt(:,nlIdx))  ).^(1/3);
        
        % rotor-averaged wind speed in the horizons
        mpc.V      = reshape(temp,Wp.turbine.N,mpc.Nh);            
         
        % Calculate cost
        uOptimize = value(U);
        dUOptimize1 = [ uOptimize(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            uOptimize(Wp.turbine.N+1:end)-uOptimize(1:end-Wp.turbine.N)];
        eOptimize1 = value(E);       
        costTestOptimize = eOptimize1'*mpc.Q*eOptimize1 + dUOptimize1'*mpc.R*dUOptimize1;
        
 %% Quadprog solver 
        tic;
        mpc = quadprog_sol(mpc,Wp,sol,xinit);
        tquadprog= toc;
        uOptimize = mpc.x;
        dUOptimize2 = [ uOptimize(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            uOptimize(Wp.turbine.N+1:end)-uOptimize(1:end-Wp.turbine.N)];
        eOptimize2 = (mpc.L_tilde*xinit + mpc.S_tilde* uOptimize - mpc.Pref(sol.k:sol.k+mpc.Nh-1));       
        costTestQuadProg = eOptimize2'*mpc.Q*eOptimize2 + dUOptimize2'*mpc.R*dUOptimize2;
 %% Lemke's Algorithm
        tic;
        mpc = sol_lemke_dU(mpc,Wp,xinit,sol);
        tlemke = toc;
        uOptimize = mpc.U+ mpc.utemp;
        dUOptimize3 = [ uOptimize(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            uOptimize(Wp.turbine.N+1:end)-uOptimize(1:end-Wp.turbine.N)];
        eOptimize3 = (mpc.L_tilde*xinit + mpc.S_tilde* uOptimize - mpc.Pref(sol.k:sol.k+mpc.Nh-1));       
        costTestLemke = eOptimize3'*mpc.Q*eOptimize3 + dUOptimize3'*mpc.R*dUOptimize3;
 %% Comparsion of the results    
        fprintf('Cost optimize: %2.3e, quadprog: %2.3e, Lemke: %2.3e \n', ...
            costTestOptimize, costTestQuadProg, costTestLemke);
        fprintf('Time optimizer: %2.3e, quadprog: %2.3e, Lemke: %2.3e, ratio: %2.3e \n', ...
            tsolver*1e3, tquadprog*1e3, tlemke*1e3,(tsolver/tquadprog));
        
    end
end
    
    %% Assign the decision variables
    Yopt          = value(Y);
    Uopt          = Yopt(mpc.MU);
    mpc.U = Uopt;
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);
    
end