function mpc = quadprog_sol(mpc,Wp,sol,xinit)

% quadprog_sol summary of this function goes here
% Solver solves for quadratic objective functions with linear constraints.
% quadprog finds a minimum for a problem specified by
% min f(x) = 0.5x'Hx+f'x
% subjected to (A)x<=b

Nh = mpc.Nh;        % prediction trajectory length
N = Wp.turbine.N;   % number turbines

% calculating the Hessian and Gradient 
if mpc.controller == 2
mpc = H_G_matrix_dU(Wp,mpc,sol,xinit);
elseif mpc.controller == 4
    mpc = H_G_Kmatrix_dU(Wp,mpc,sol,xinit);
end

% Subjected to  (Aco_)x<=bco_
Ac0_ = [eye(Nh*N); -eye(Nh*N)];
bc0_ = [mpc.Ulim_lower; - mpc.Ulim_upper];
bc   = -bc0_ + Ac0_ * mpc.utemp;

% Turning into QP and solve it%%%%%%%%%
myoptions = optimset('Display','off'); 
[uval ,fval] = quadprog(mpc.H,mpc.g,-Ac0_,(-bc0_ + Ac0_ * mpc.utemp),[],[],[],[],[],myoptions );

%[uval ,fval]= quadprog(mpc.H,mpc.g,-Ac0_,bc);%myquadprog(mpc.H,mpc.g,-Ac0_,bc);

mpc.x = uval + mpc.utemp;
mpc.uval = uval;
mpc.fval = fval;
