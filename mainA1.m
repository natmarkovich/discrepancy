%% Calculation of intervals estimate with "plateau-finding" algorithm A1 by M. Ferreira (2018)
%% Heuristic tools for the estimation of the extremal index: a comparison of methods. 
%% REVSTAT - Statistical Journal,  16(1), 115-136.
%% Paper Markovich and Rodionov 2018
clc;
clear;

%% Simulated stochastic processes:
% 1 Moving Maxima process, 2 uniform AR(1) , 3 ARMAX process 4 Cauchy AR(1)
% 5 AR(2) with Pareto noise % 6 Pareto MA(2) % 7 GARCH(1,1)
% AR(1) with negative correlation
%distribution_type = 1;
%distribution_type = 2;
%distribution_type = 3;
%distribution_type = 4;
%distribution_type = 5;
%distribution_type = 6;
distribution_type = 7;
%distribution_type = 9;

sample = 1000; % sample of repeating
n = 5000; % number of elements

%Parameters of the MM process 1
MM_parametr.m = 3;
 MM_parametr.a1 = 0.5;
 MM_parametr.a2 = 0.3;
 MM_parametr.a3 = 0.15;
 MM_parametr.a4 = 0.05;
 MM_parametr.theta = 0.5; % Extremal Index
MM_parametr.a1 = 0.8;
MM_parametr.a2 = 0.12;
MM_parametr.a3 = 0.05;
MM_parametr.a4 = 0.03;
MM_parametr.theta = 0.8; % Extremal Index

%% Parameters of the AR(1) with uniform noise
%AR_parametr.r = 2;
%AR_parametr.theta = 0.5;% Extremal Index
AR_parametr.r = 5;
%AR_parametr.theta = 0.8;% Extremal Index
AR_parametr.theta = 1-1/AR_parametr.r.^2;% Extremal Index
%% Parameters of the AR(1) with Cauchy noise
ARC_parametr.theta = 0.3;% Extremal Index
%% Parameters of the AR(2) with Pareto noise
AR2_parametr.alpha = 2; % Tail index
AR2_parametr.theta = 0.25;% Extremal Index
%% Parameters of the MA(2) with Pareto noise
MA2_parametr.alpha = 2; % Tail index
MA2_parametr.p = 1/2.^0.5; % Parameter p
MA2_parametr.q = 1/2.^0.5; % Parameter q
MA2_parametr.theta = 0.5;% Extremal Index
MA2_parametr.p = 1/3.^0.5; % Parameter p
MA2_parametr.q = 1/6.^0.5; % Parameter q
MA2_parametr.theta = 2/3;% Extremal Index
%% Parameters of the ARMAX
ARMAX_parametr.alfa = 0.75;
ARMAX_parametr.theta = 0.25;% Extremal Index
%ARMAX_parametr.alfa = 0.25;
%ARMAX_parametr.theta = 0.75;% Extremal Index
%ARMAX_parametr.alfa = 0.5;
%ARMAX_parametr.theta = 0.5;% Extremal Index

%% Parameters of the GARCH(1,1) with gaussian noise
GARCH_parametr.alpha = 1/(10.^6);
GARCH_parametr.lambda = 0.25;
GARCH_parametr.beta = 0.7;
GARCH_parametr.theta = 0.447;% Extremal Index

%% Parameters of the Lindley for M/M/1 with lambda<mu
LINDLEY_parametr.lambda = 0.5;
LINDLEY_parametr.mu = 0.7;
LINDLEY_parametr.theta = 0.082;

% X includes all samples we generate
X = generation_X(distribution_type, n, sample, MM_parametr,...
AR_parametr, ARMAX_parametr, AR2_parametr,MA2_parametr,GARCH_parametr,LINDLEY_parametr);

 s_1=0; s_2=0; s_3=0; n_1=0; j_1=1; th_1=0;
 s_1
for j = 1 : sample
    j
       [theta_1,  flag_s] = intervalsestimatorA1(X(:, j));
              flag_s
       if(flag_s ==0)
           th_1(j_1)=theta_1;  
           j_1=j_1+1;
               s_1=s_1+theta_1;
               s_1
               
                      else
       n_1=n_1+1;
                    end
     
end 
%th_1
% Standard Error=Standard deviation/sqrt(n)
% Bias=mean_bias over resamples - theta
% Root MSE
%ARMAX_parametr.theta=AR_parametr.theta;
%ARMAX_parametr.theta=ARC_parametr.theta;
%ARMAX_parametr.theta=AR2_parametr.theta;
%ARMAX_parametr.theta=MA2_parametr.theta;
%ARMAX_parametr.theta=MM_parametr.theta;
ARMAX_parametr.theta=GARCH_parametr.theta;
ARMAX_parametr.theta
m=sample-n_1;
m
%s_1/(sample-n_1)
 bias_1=0; 
bias_1=s_1/m-ARMAX_parametr.theta;
bias_1
%var(th_1)
rmse_1=(var(th_1)+bias_1.^2).^0.5;
rmse_1 % root mse


