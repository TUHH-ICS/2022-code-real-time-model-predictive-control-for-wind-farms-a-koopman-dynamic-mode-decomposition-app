clc; clear; close all;

%% 1. Create or load open loop test data with WFSim
addpath('WFSimCode');

% Inputs
R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = 0; %refstairs: Set reference to stairs
measured = 1; %measured (for Koopman model): Use measured values as feedback
KoopmanStates =  24; %KoopmanStates (for Koopman model):  Number of Koopman states
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24
controller = 0; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
ControlSetStr = 'sowfa_2turb_alm_turbl_AllComb';
Vinf = 8;
vinfstr = '';

% Run WFSim_demo
WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);

%% 2. Generate Koopman Sys Id for windfarm controller
addpath('KoopmanIODMD')

% Inputs
yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control
filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb.mat';
filenameVal = 'Vinf8dot5_sowfa_2turb_alm_turbl.mat';
noStates = [6,12,12,14,16,18,24]; %number of Koopman states
useVal = 0; % use extra data for validation
percentTrain = .6;

% Run Koopman main function for WFSim simulation environent
MainWfSim(yawmode,filenameId,filenameVal,noStates,useVal,percentTrain);


%% 3. Evaluate quality of Koopman Sys ID in WFsim in closed loop with MPC
% Inputs

R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = 0; %refstairs: Set reference to stairs
measured = 1; %measured (for Koopman model): Use measured values as feedback
controller = 2; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl';
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24
KoopmanStates = 12; %KoopmanStates (for Koopman model):  Number of Koopman states
% Run WFSim_demo

WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf); % first for measured values
measured = 0;
%KVec = [2,4,6,10,12,12,14,16,18,24];
PolyVec = zeros(size(noStates)); % 1: use only polynomial, 0: otherwise
idxPolyOn = find(noStates == 12, 1, 'first');
PolyVec(idxPolyOn) = 1;
for idx = 1: length(noStates) % first for estimated values
    KoopmanStates = noStates(idx);
    PolyLiftingFunction = PolyVec(idx);
    try
        [~,JR(idx),JQ(idx)] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf);
    catch me
        save(['Error',num2str(idx),'.mat'],'me');
    end
    
end
idx = 2;
KoopmanStates = noStates(idx);
PolyLiftingFunction = PolyVec(idx);
[~,JR(idx),JQ(idx)] = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf);
