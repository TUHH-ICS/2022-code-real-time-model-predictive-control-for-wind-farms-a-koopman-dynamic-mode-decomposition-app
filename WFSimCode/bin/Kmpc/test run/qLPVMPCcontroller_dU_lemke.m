function [CT_prime,phi,mpc] = qLPVMPCcontroller_dU_lemke(sol,Wp,mpc)

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

        % build wind farm model
        mpc  = wfmodel(sol,Wp,mpc,nlIdx);
        
        % build matrices horizon
        mpc  = Matsys(Wp,mpc,mpc.Nh);
         
        % solve the optimization problem and calculate time
        tic
        mpc = sol_lemke_dU(mpc,Wp,xinit,sol);
        mpc.displaytime = toc;
        
        % state space model of the wind farm
        aValue_approxLemke = mpc.U + mpc.utemp;
        X    = mpc.AA*xinit + mpc.BBt*aValue_approxLemke;   
        Y    = mpc.CC*X;                                    %yopt
        P    = reshape(Y(mpc.MP),Wp.turbine.N,mpc.Nh);      %Power output
        
        % tracking error value
        Etemp_L    = mpc.Pref(sol.k:sol.k+mpc.Nh-1) - sum(P)';
        
        % change in input error value
        d_utemp2   = [aValue_approxLemke(1:Wp.turbine.N)-sol.turbine.CT_prime ; ...
            aValue_approxLemke(Wp.turbine.N+1:end) - aValue_approxLemke(1:end-Wp.turbine.N)];
        
        % Calculate the cost function
        mpc.costValue_Lemke = Etemp_L'*mpc.Q*Etemp_L + d_utemp2'*mpc.R*d_utemp2;

    end
end  
    %% Assign the decision variables
    Yopt          = Y;
    Uopt          = Yopt(mpc.MU);
    temp          = reshape(Uopt,[Wp.turbine.N,mpc.Nh]);
    
    CT_prime      = temp(:,1);              % first action horizon
    phi           = zeros(Wp.turbine.N,1);  % Yaw input is neglected
    

end
