function [Inputs, Outputs, Deterministic,scalingfactors] = preprocessdmdid(beg,dataSOFWA,pitchmode,rho,meanY1,meanY2,meanX1,meanX2)


%% Mean not subtracted for absolute power
if nargin <= 4
   meanY1 = 5.352*10^6; %yaw -10
   meanY2 = 0.9512*10^6; %yaw -10
   meanX1 = 7.14; %pitch 0
   meanX2 = 4.027; %yaw -10
end
%% EVALUATE RELEVANT STATES TO BE UED
idx = 1;
time = dataSOFWA{idx}.time ;
pitch = dataSOFWA{idx}.pitch;
powerGenerator = dataSOFWA{idx}.powerGenerator;
rotSpeed = dataSOFWA{idx}.rotSpeed;
rotorAzimuth = dataSOFWA{idx}.rotorAzimuth;
nacelleYaw = dataSOFWA{idx}.nacelleYaw;



X1 = (rotSpeed(end-beg*10:1:end,1)-meanX1).^2;
X2=(rotSpeed(end-beg*10:1:end,2)-meanX2).^2;

%rotor speed of turbine 2 as second output
% X1=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10);
% X2=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10);
% X3=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10).^2;
% X4=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10).^2;

[X1] = resampleedgeeffect(X1,10);
[X2] = resampleedgeeffect(X2,10);
%[X3] = resampleedgeeffect(X3,10);
%[X4] = resampleedgeeffect(X4,10);

%s4=var(X1);
%s5=var(X2);
%X1=X1./var(X1);
%X2=X2./var(X2);

%%
% X3=X1.^2;
% X4=X2.^2;
%s6=var(X3); %scaling factor of first deterministc states
%s7=var(X4);
%X3=X3./var(X3);
%X4=X4./var(X4);

% X3=resample(rotSpeed(end-beg*10:1:end,1),1,10).^2;
% X4=resample(rotSpeed(end-beg*10:1:end,2),1,10).^2;
% X3=X3./var(X3);
% X4=X4./var(X4);

%% INPUTS:
if pitchmode==0
    
    steadyyaw=260; %yaw angle steady state (offset)
    %U1=detrend(nacelleYaw(end-beg*10:1:end,1)');
    U1 = nacelleYaw(end-beg*10:1:end,1)';
    U1 = U1-steadyyaw;
    U1 = resampleedgeeffect(U1,10);
    %s3=var(U1);
    %U1=U1./var(U1); %scaling
    
    Inputs = U1;
    
    % transformation to make input more linear
    % Inputs=(1-cosd(Inputs))*6;
    
    %scalingfactors=[s1;s2;s3;s4;s5;s6;s7];
    scalingfactors=[1;1;1;1;1;1;1];
    %scalingfactors=0;
    
elseif pitchmode==1
    
    %% MBC: Multi-Blade Coordinate transformation
    %A directional thrust force  can be accomplished by implementinf MBC
    %transformation, and decoupling/proejcting the blade loads in a non
    %-rotating reference frame
    
    % As a result, the measured out-of plane blade root bending moments M(t)
    % --> [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)] are
    % projected onto a non rotating reference frame --> PITCH
    
    Nturb=2;
    Offset=-8.4*2;
    
    for index=1:1:length(time)
        %for each time instant INDEX get (for each turbine ij below) the 3
        %out-of-plane blade root bending moments, corresponding to the
        %three columns given a certain line
        
        for ij=1:1:Nturb
            
            Azimuth=rotorAzimuth(index,ij);
            
            PITCH=([1/3 1/3 1/3;
                2/3*cosd(Azimuth+Offset) 2/3*cosd(Azimuth+120+Offset) 2/3*cosd(Azimuth+240+Offset);
                2/3*sind(Azimuth+Offset) 2/3*sind(Azimuth+120+Offset) 2/3*sind(Azimuth+240+Offset);])*...
                [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)];
            
            Pitch1(ij)=PITCH(1);
            Pitch2(ij)=PITCH(2);
            Pitch3(ij)=PITCH(3);
            
            %3 Matrixes containing the different bending moments where each
            %line has the turbine number and each column the time instant
            %INDEX
            PPitch1(ij,index)=PITCH(1);
            PPitch2(ij,index)=PITCH(2);
            PPitch3(ij,index)=PITCH(3);
            
        end
    end
    
    %U1=resample(detrend(PPitch2(1,end-750*10:1:end)),1,10);
    %s3=var(U1);
    %U1=U1./var(U1);
    %U2=resample(detrend(PPitch3(1,end-750*10:1:end)),1,10);
    %s4=var(U2);
    %U2=U2./var(U2);
    U1 = PPitch1(1,end-beg*10:1:end);
    U1 = resampleedgeeffect(U1,10);
    Inputs= [U1];
    
    scalingfactors=[1;1;1;1;1;1;1];
    
end


%% OUTPUTS
%the rotor speeds of the two turbines are defined as outputs of the wind turbine system

%OUTPUT POWER
% Y1=resample(detrend(powerGenerator(end-750*10:1:end,1)*1e-6'),1,10);
% Y2=resample(detrend(powerGenerator(end-750*10:1:end,2)*1e-6'),1,10);

% Y1=detrend(powerGenerator(end-beg*10:1:end,1)');
% Y2=detrend(powerGenerator(end-beg*10:1:end,2)');

% meanY1=5.485*10^6;%pitch 0
Y1= powerGenerator(end-beg*10:1:end,1)/rho  -meanY1;

% meanY2=0.7728*10^6; %pitch 0
Y2=powerGenerator(end-beg*10:1:end,2)/rho - meanY2;

[Y1] = resampleedgeeffect(Y1*10^-6,10); %rotor speed of turbine 1 as first output
[Y2] = resampleedgeeffect(Y2*10^-6,10);

% TRANSFORMATION
%  Y1=Y1.*Inputs.^(1.88);

%s1=1;
%s2=1;
%s1=var(Y1);
%s2=var(Y2);

%Y1=Y1./var(Y1);
%Y2=Y2./var(Y2);

%OUTPUT ROTOR SPEED
% Y1=X1;
% Y2=X2;

Outputs=[Y1;Y2];



%% DETERMINISTIC STATES

Deterministic=[X1; X2];
%Deterministic={};

%% Storing of scaling factors

% the scaling factors are stored to be used for validation later, to ensure
% that identification and validation data are in the same order of
% magnitude

%scaling factor for output: S1
%scaling factor for output: S2

%sacling factor for input: S3

%scaling factor for deterministic states: S4
%scaling factor of second deterministc state: S5






