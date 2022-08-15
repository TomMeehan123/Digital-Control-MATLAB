%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename: Assignment 3 - Assignment_3.m
%
% purpose: Use MATLAB to solve for state space equations
%
% created by: Tom Meehan
% created on: 25 May 2022
%
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------
clc
clear
clear all
close all

%% Question 1

% Create the state space matricies in the continious time domain (A,B,C,D)

A = [0 1; 0 -3.75];
B = [0; 1];
C = [1 0];
D = [0];

% Create a state space model using the matricies

Gss = ss(A,B,C,D);

%% Question 2

% Create the state space equations in the discrete time domain
T = 1;

[phi,gama] = c2d(A,B,T);

Gzs = ss(phi,gama,C,D,T);

[num,den] = ss2tf(phi,gama,C,D);

Gz = tf(num,den,T);

%% Question 5

n = 21; % samples
k = 0:n-1;
u = [1 zeros(1,n-1)]; % input u(k) = impulse
 
% Determine the impulse response for the discrete time system
num_GZ = [0.1972 0.06317 0]; % numerator for Gz
den_GZ = [1 -1.02352 0.02352]; % denominator for Gz
yk = filter(num_GZ,den_GZ,u); % impulse using the filter command
 
% Determine the step response for the discrete time system
num_GZ_Step = [0.1972 0.06317 0 0]; % numerator for Gz
den_GZ_Step = [1 -2.02352 1.04704 -0.02352]; % denominator for Gz
yk_step = filter(num_GZ_Step,den_GZ_Step,u); % impulse using the filter command

figure(1)
subplot(2,1,1)
hold on
stairs(k,yk,'b') % stairs function plots the impulse
title('Impulse Response for the Discrete Time System')
xlabel('Sample (y(k))')
ylabel('Amplitude')

figure(2)
subplot(2,1,1)
hold on
stairs(k,yk_step,'b') % stairs function plots the impulse
title('Step Response for the Discrete Time System')
xlabel('Sample (y(k))')
ylabel('Amplitude')
 
% Use impulse command to verify results
figure(1)
subplot(2,1,2)
hold on
impulse(Gss,'b')
title('Impulse Response for the Continious Time System')

figure(2)
subplot(2,1,2)
hold on
step(Gss,'b')
title('Step Response for the Continious Time System')
