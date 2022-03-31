function mpc = MPCinitKoopman(mpc,sol,Wp,NN)

% controller models
%KWFModel = load('Method-1_Koopman1_states4_stateNameUr1;Ur2;Ur1.^2;Ur2.^2.mat');
KWFModel = load('Method-1_Koopman1_states8_stateNameUr1;Ur2;Ur1.^2;Ur2.^2;Ur1.^3;Ur2.^3;dUr1.^3;dUr2.^3.mat');
mpc.A = KWFModel.sys_red{1,1}.A;   
mpc.B = KWFModel.sys_red{1,1}.B;
mpc.C = KWFModel.sys_red{1,1}.C(3,:)+KWFModel.sys_red{1,1}.C(4,:); % to have only P1 and P2 as output
mpc.D = KWFModel.sys_red{1,1}.D(3,:)+KWFModel.sys_red{1,1}.D(4,:); %outputs =[Ur1 Ur2 PT1 PT2 FT1 FT2]
mpc.nx = size(mpc.A,1);
mpc.Xinit = KWFModel.xo;
