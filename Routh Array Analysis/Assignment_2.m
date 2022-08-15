%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename: Assignment 2 - Assignment_2.m
%
% purpose: Use MATLAB to solve for a ZoH system
%
% created by: Tom Meehan
% created on: 25 March 2022
%
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------
clc
clear
clear all
close all

%% Question 1
% Determine the ZOH
T = 0.1; % sampling period, T
Gs = tf(1.53125,[1 3.0625 0]); % G(s)
Gz = c2d(Gs,T); % ZOH equivalent, G(z)

% Verify the partial fractions 

num = [1];
den = [1 3.0625 0 0];

[r,p,k] = residue(num,den); % Calculate the residues A,B,C

Ku = 42.27;
Kz = Ku/8;
GclZ = feedback(Gz*Kz,1);

figure(1)
pzmap(GclZ)

%% Question 3
% plot the root locus for Gz which is the start point when Kp = 0
figure(2)
rlocus(Gz)

%% Question 4
n = 51; % samples
k = 0:n-1;
u = [1 zeros(1,n-1)]; % input u(k) = impulse

% Determine the impulse response
num_GZ = [0.03646 0.03297 0]; % numerator for Gz
den_GZ = [1 -1.6997 0.7692]; % denominator for Gz
yk = filter(num_GZ,den_GZ,u); % impulse using the filter command

figure(3)
hold on
stairs(k,yk*10,'b') % stairs function plots the impulse

% Use impulse command to verify results
figure(4)
hold on
impulse(GclZ)

% Determine the step response
hold on
num_Sz = [0 0.03646 0.03297 0]; % numerator for Gz step function
den_Sz = [1 -2.6997 2.4689 -0.7692]; % denominator for Gz step function

sk = filter(num_Sz,den_Sz,u); % step function for Gz
figure(3)
hold on
stairs(k,sk,'r')
xlabel('Sample Period (k)')
ylabel("Amplitude")
title("Impulse Response & Step Response")
legend('Impulse Response','Step Response')
axis([0 50 -0.4 1.6])

% Use step command to verify results
figure(4)
hold on
step(GclZ)
title("Impulse Response & Step Response")
legend('Impulse Response','Step Response')

%% Question 5
num_GZ_ramp_Input = [0 T 0]; % numerator for Gz
den_GZ_ramp_Input = [1 -2 1]; % denominator for Gz
yk_ramp_In = filter(num_GZ_ramp_Input,den_GZ_ramp_Input,u); % impulse using the filter command

figure(5)
hold on
stairs(k,yk_ramp_In,'b') % stairs function plots the impulse

num_GZ_ramp_Response = [0 0 0.03646*T 0.03297*T 0]; % numerator for Gz
den_GZ_ramp_Response = [1 -3.6992 5.1676 -3.2376 0.7692]; % denominator for Gz
yk_ramp = filter(num_GZ_ramp_Response,den_GZ_ramp_Response,u); % impulse using the filter command

figure(5)
hold on
stairs(k,yk_ramp,'r') % stairs function plots the impulse
title("Ramp Input & Ramp Response")
legend('Ramp Input','Ramp Response')