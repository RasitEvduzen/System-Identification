%   13-Apr-2018 14:16:26
%   Discrete Time linear kalman filter
%   Rasit EVDÜZEN
%% Discrete linear kalman filter Design
clc,clear all,close all,warning off;

% Sampling Time and Period
Ts = 1e-2;
T  = 1e1;

% System initials
x = [-1;1];     % System initial state
N = length(x);  % Number of state
A = [0 1; -2 -2];   
B = [0;1];
C = [1 0];
D = 0;

% Discretization of system
sys = ss(A,B,C,D);
s1 = c2d(sys,Ts,'Tustin');
Ad = s1.a; Bd = s1.b; Cd = s1.c;
% Ad = [0 1; 0.5 0.2]; Bd = (exp(A*Ts)-eye(N))*inv(A)*B; Cd = C;

% Kalman initials

xkf = [0.1;0.1];   % Kalman Filter initial State (Tuning Variable)
Pkf = 1*eye(N); 
Q = 0.001*[1 1; 1 1];       % Process noise (State Noise)
R = 0.5;                    % Measurement noise (Output Noise)

% Main Loop
for n = 1:T/Ts
    % System integration
    t(n) = n*Ts;
    u(n) = 5*cos(1*t(n));
    Noise(:,n) = [sqrt(Q(1,1))*randn;sqrt(Q(2,2))*randn];  % Processes Noise
    x(:,n+1) = Ad * x(:,n) + Bd*u(n) + Noise(:,n);  % Euler integration
    ynoise(n) = sqrt(R)*randn;
    y(n) = Cd*x(:,n)+ynoise(n);   % Measurement Noise
    
    % Time Update
    xp(:,n+1) = Ad*xkf(:,n)+Bd*u(n);
    z(n) = Cd*xp(:,n+1);
    Pp = Ad*Pkf*Ad'+Q;
    
    % Posteriori Update
    K = Pp*Cd'*pinv(Cd*Pp*Cd'+R);   % Kalman Gain Calculate
    xkf(:,n+1) = xp(:,n+1) + K*(y(n)-z(n));
    Pkf = (eye(N)-K*Cd)*Pp;
    
    % Estimation Error
    e(:,n) = x(:,n) - xkf(:,n);
        
end
% SNR Calculation
Ps = y*y'*Ts/T;
Pn = ynoise*ynoise'*Ts/T;
SNR = 10*log(Ps/Pn);

Ps2 = x(2,:)*x(2,:)'*Ts/T;
Pn2 = Noise(2,:)*Noise(2,:)'*Ts/T;
SNR2 = 10*log(Ps2/Pn2);

% RMS Calculation
RMSE_x2 = e(2,:)*e(2,:)'*Ts/T;



%% Plot Data
subplot(221)
plot(t,x(1,2:end));
hold on,grid minor
plot(t,xkf(1,2:end),'r--');
legend('Real State x1','Estimate State x1')

subplot(222)
plot(t,x(2,2:end));
hold on,grid minor
plot(t,xkf(2,2:end),'r--');
legend('Real State x2','Estimate State x2')

subplot(223)
plot(t,y)
hold on,grid minor
plot(t,z,'r--')
legend('System Output','Kalman Output')






