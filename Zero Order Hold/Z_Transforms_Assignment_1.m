%--------------------------------------------------------------------------
% University of Limerick - Dept. of Electronic and Computer Engineering
%--------------------------------------------------------------------------
% filename: Assignment 1 - Z_Transforms_Assignment_1.m
%
% purpose: Use MATLAB to solve for a ZoH system
%
% created by: Tom Meehan
% created on: 24 Feburary 2022
%
%--------------------------------------------------------------------------
% Copyright 2021 University of Limerick
%--------------------------------------------------------------------------
clc
clear
clear all
close all

%% Question 1
% Create a variable for s
syms s

% Verify the partial fractions 

num = [1.75];
den = [1 1.75 0 0];

[r,p,k] = residue(num,den) % Calculate the residues A,B,C

x = partfrac(1/(s^3+1.75*s^2)); % Calculate the partial fraction expansion
pretty(x) % Pretty print

% Determine the ZOH
T = 0.1; % sampling period, T
Gs = tf(1.75,[1 1.75 0]); % G(s)
Gz = c2d(Gs,T) % ZOH equivalent, G(z)

%% Question 2
figure(1)
pzmap(Gz) % Pole Zero Plot of Gz

Zg = zero(Gz) % Calculate the zeros for GZ
Pg = pole(Gz) % Calculate the poles for GZ

%% Question 3
KZ = 1.75; % Assign value for gain Kz

GclZ = feedback(Gz*KZ,1) % Calculate the closed loop transfer function
figure(2)
pzmap(GclZ) % Pole Zero Plot of GclZ

Zgcl = zero(GclZ) % Calculate the zeros for GclZ
Pgcl = pole(GclZ) % Calculate the poles for GclZ

Sgcl = abs(Pgcl) % Determine if GclZ is stable
%% Question 5
n = 120; % samples
k = 0:n-1;
u = [1 zeros(1,n-1)]; % input u(k) = impulse

num_GZ = [0 0.008 0.00768]; % numerator for Gz
den_GZ = [1 -1.839 0.8395]; % denominator for Gz
yk = filter(num_GZ,den_GZ,u); % impulse using the filter command

figure(3)
hold on
stairs(k,yk,'b') % stairs function plots the impulse

num_GclZ = [0 0.014 0.01344]; % numerator for Gclz
den_GclZ = [1 -1.825 0.8531]; % denominator for Gclz
yk_cl = filter(num_GclZ,den_GclZ,u); % impulse using the filter comman

stairs(k,yk_cl,'r') % stairs function plots the impulse
xlabel('Sample Period (k)')
ylabel("Amplitude")
title("Impulse Response")
legend('G(Z) Plant Model','Gcl(Z) Closed Loop Model')

figure(4)
hold on
num_Sz = [0 0.008 0.00768 0]; % numerator for Gz step function
den_Sz = [1 -2.839 2.6785 -0.8395]; % denominator for Gz step function

sk = filter(num_Sz,den_Sz,u); % step function for Gz
stairs(k,sk,'b')

num_SclZ = [0 0.014 0.01344 0]; % numerator for Gclz step function
den_SclZ = [1 -2.825 2.6781 -0.8531]; % denominator for Gclz step function

sk_cl = filter(num_SclZ,den_SclZ,u); % step function for Gclz
stairs(k,sk_cl,'r')
xlabel('Sample Period (k)')
ylabel("Amplitude")
title("Step Response")
legend('G(Z) Plant Model','Gcl(Z) Closed Loop Model')