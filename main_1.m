%% Discrepancy method based on a normalized von Mises - Smirnov statistic W2 built by k largest order statistics of the sample 
%% Source: Algorithm 1 in Markovich and Rodionov 2022 Threshold selection for extremal index estimation. arXiv:2009.02318v1
clc;
clear;

%% Simulated stochastic processes:
% 1 Moving Maxima process MM, 2 uniform AR(1) process, 3 ARMAX process 4
% Cauchy AR(1) process
% 5 AR(2) process with Pareto noise % 6 Pareto MA(2) process % 7 GARCH(1,1) process %9 uniform
% AR(1) process with negative correlation % 10 Lindley process
%
% distribution_type shows the number of selected stochastic process
%distribution_type = 1;
%distribution_type = 2;
distribution_type = 3;
%distribution_type = 4;
%distribution_type = 5;
%distribution_type = 6;
%distribution_type = 7;
%distribution_type = 9;
%distribution_type = 10;

sample = 5000; % sample of repeating
n = 5000; % number of elements

%Parameters of the MM process 1
MM_parametr.m = 3;
 MM_parametr.a1 = 0.5;
 MM_parametr.a2 = 0.3;
 MM_parametr.a3 = 0.15;
 MM_parametr.a4 = 0.05;
 MM_parametr.theta = 0.5; % Extremal Index
 %Parameters of the MM process 2
%MM_parametr.a1 = 0.8;
%MM_parametr.a2 = 0.12;
%MM_parametr.a3 = 0.05;
%MM_parametr.a4 = 0.03;
%MM_parametr.theta = 0.8; % Extremal Index

%% Parameters of the AR(1) process with uniform noise 1
%AR_parametr.r = 2;
%AR_parametr.theta = 0.5;% Extremal Index
%% Parameters of the AR(1) process with uniform noise 2
AR_parametr.r = 5;
%AR_parametr.theta = 0.8;% Extremal Index
%% Parameters of the AR(1) process with negative correlation
AR_parametr.theta = 1-1/AR_parametr.r.^2;% Extremal Index
%% Parameters of the AR(1) with Cauchy noise
ARC_parametr.theta = 0.3;% Extremal Index
%% Parameters of the  AR(2) with Pareto noise
AR2_parametr.alpha = 2; % Tail index
AR2_parametr.theta = 0.25;% Extremal Index
%% Parameters of the MA(2) with Pareto noise 
MA2_parametr.alpha = 2; % Tail index
%MA2_parametr.p = 1/2.^0.5; % Parameter p
%MA2_parametr.q = 1/2.^0.5; % Parameter q
%MA2_parametr.theta = 0.5;% Extremal Index
MA2_parametr.p = 1/3.^0.5; % Parameter p
MA2_parametr.q = 1/6.^0.5; % Parameter q
MA2_parametr.theta = 2/3;% Extremal Index
%% Parameters of the ARMAX process
%ARMAX_parametr.alfa = 0.75;
%ARMAX_parametr.theta = 0.25;% Extremal Index
%ARMAX_parametr.alfa = 0.25;
%ARMAX_parametr.theta = 0.75;% Extremal Index
ARMAX_parametr.alfa = 0.8;
ARMAX_parametr.theta = 0.2;% Extremal Index
%ARMAX_parametr.alfa = 0.5;
%ARMAX_parametr.theta = 0.5;% Extremal Index

%%Parameters of the GARCH(1,1) with gaussian noise
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
 s_1=0; s_2=0; s_3=0; s_4=0; n_1=0; j_1=1; th_1=0;
 s_1
for j = 1 : sample
    j
       [theta_1, theta_2, theta_3, flag_s] = intervalsestimator(X(:, j)); 
              flag_s
       if(flag_s ==0)
           th_1(j_1)=theta_1; th_2(j_1)=theta_2; th_3(j_1)=theta_3;  
           j_1=j_1+1;
               s_1=s_1+theta_1;
               s_1
               s_2=s_2+theta_2;
       s_2
       s_3=s_3+theta_3;
       s_3
              else
       n_1=n_1+1;
                    end
        %   [thr_val_q, theta, w_plus] = intervalsestimator(X(:, j));
 %   for u_c=1:length(thr_val_q)
 %       u=prctile(X(:, j),thr_val_q(u_c));
 %       u
 %       [theta, Y]=intervalsestimatoru(X(:, j),u);
                      % Calculation w2 by function count_w
       
 %       W_1 = @(u_c) count_W(Y, theta(u_c));
 %       u_W_1 =@(u_c) (W_1(u_c) - 1.49).^2;
 %   U_1 = fminbnd(u_W_1, 0, 500);
 %   W_1(U_1);
 %   end
    
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
%ARMAX_parametr.theta=GARCH_parametr.theta;
%ARMAX_parametr.theta=LINDLEY_parametr.theta;
ARMAX_parametr.theta
m=sample-n_1;
m
%s_1/(sample-n_1)
 bias_1=0; bias_2=0;bias_3=0;
bias_1=s_1/m-ARMAX_parametr.theta; % Bias
bias_1
%var(th_1)
std_1=std(th_1)/m.^0.5;
std_1
rmse_1=(var(th_1)+bias_1.^2).^0.5;
rmse_1 %Root mean squared error (MSE)
bias_2=s_2/m-ARMAX_parametr.theta;
bias_2
%var(th_2)
std_2=std(th_2)/m.^0.5; % Standard deviation
std_2
rmse_2=(var(th_2)+bias_2.^2).^0.5;
rmse_2
bias_3=s_3/m-ARMAX_parametr.theta;
bias_3
std_3=std(th_3)/m.^0.5;
std_3
rmse_3=(var(th_3)+bias_3.^2).^0.5;
rmse_3
%v_1=var(th_1);
%v_1

