function turbInputSet = controlSet_sowfa_2turb_yaw_alm_combined2(Wp)

if ~nargin % enable running this as a script
    Wp.turbine.Crx = nan(2,1);
end

turbInputSet1 = controlSet_sowfa_2turb_yaw_alm_turbl;
turbInputSet2 = controlSet_sowfa_2turb_yaw_alm_uniform;

k1 = 875; %te = turbInputSet1.t(k1);
k2 = 60; %t0 = turbInputSet1.t(k2);

% turbInputSet.t = [turbInputSet1.t(1:k1), turbInputSet2.t(k2:end-1) +(te +1 -t0)];
tempPhi = [turbInputSet1.phi(:,1:k1), turbInputSet2.phi(:,k2:end)];
turbInputSet.t = 0: length(tempPhi)-1;

temp = [turbInputSet1.CT_prime(:,1:k1), ...
    turbInputSet2.CT_prime(:,k2:end)];

tfD = c2d(tf(1,[20,1]),1);
B = tfD.Numerator{:};
A = tfD.Denominator{:};

rng(0);
% X = randn(1,k1);
% X1 = randn(1,k1);
% Xfilt(1,:) = filter(B,A,X);
% Xfilt(2,:) = filter(B,A,X1);
% rnd1 = [Xfilt, zeros(2,length(turbInputSet2.CT_prime(:,k2:end)))];

X = randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end)))*2;
X1 = randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end)))*2;
Xfilt(1,:) = filter(B,A,X);
Xfilt(2,:) = filter(B,A,X1);
rnd1 = Xfilt;


tfD1 = c2d(tf(1,[30,1]),1);
B1 = tfD1.Numerator{:};
A1 = tfD1.Denominator{:};


X = (randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end))))*30;
X1 = (randn(1,k1+length(turbInputSet2.CT_prime(:,k2:end))))*30;
Xfilt1(1,:) = filter(B1,A1,X);
Xfilt1(2,:) = filter(B1,A1,X1)*0;
rnd2 = Xfilt1;

% stdT = std(temp);
turbInputSet.CT_prime  = min(max((temp + rnd1),0.2),2) ;

turbInputSet.phi =  tempPhi + rnd2;
turbInputSet.interpMethod = 'lin'; % Linear interpolation over time

if length(Wp.turbine.Crx) ~= size(turbInputSet.phi,1)
    error('Number of turbines in layout does not match your controlSet.');
end
end