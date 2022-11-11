function [theta_1, theta_2, theta_3, flag_s] = intervalsestimator(X)
%Intervals estimator with discrepancy method as a threshold selection tool.
%The intervals estimator is used as a pilot estimator of
%extremal index to compute k largest order statistics (see, Algorithm 1
%(step (2)) in Markovich and Rodionov (2022) arXiv:2009.02318v1
%Intervals Estimator estimator==============================================
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%IE estimator==============================================
X_l= length(X); % Length of the sample we need for Intervals Estimator
thr_val_q = [90 90.5 91 91.5 92 92.5 93 93.5 94 94.5 95 95.5 96 96.5 97 97.5 98 98.5 99 99.5]; % percentage quantiles
%thr_val_q = [60 60.5 61 61.5 62 62.5 63 63.5 64 64.5 65 65.5 66 66.5 67 67.5 68 68.5 69 69.5]; % percentage quantiles
L=0; % The number of inter-cluster times
%s=0.05; % proportion of largest inter-cluster times
s=0.56;
%s=0.3;
%s=0.8;
%s=0.05;
eps=0.01; %  Upper deviation between the statistic w2 and sigma_1 (see, formula (13) in Markovich and Rodionov (2022))
%eps=0.02;
sigma_1=0.05; % The mode of the statistic w2 (see, formula (14) in Markovich and Rodionov (2022))
sigma_1=1.49; % The 99.98% quantile of the statistic w2 used in formula (16) as sigma_2
f_m=0;        % Number of unseccesful solutions of the discrepancy equation
flag_s=0; % The indicator of sample when  there is no solution of the discrepancy equation (14) in Markovich and Rodionov (2022)
theta_1=0; theta_2=0; theta_3=0; %Initial values of statistics in formula (15) by Markovich and Rodionov (2022)
for u_c=1:length(thr_val_q)
        u=prctile(X,thr_val_q(u_c)); % The threshold selection as quantiles of levels thr_val_q (see, Algorithm 1, step 1)
        flag=0;
%% Computation of inter-exceedance times T_i(u) (see, Algorithm 1, step 1)       
    T_t=0;T=0;Y=0;
    
    n=0; 
    
    for i=1:X_l
        
        if (flag==0)
            if (X(i)>u)
                flag=1;
                T_t=i;
            end
            
        elseif (flag==1)
            if (X(i)>u)
                n=n+1;
                T(n)=i-T_t;
                T_t=i;
            end
        end
    end
    n
    % n is the number of exceedances over threshold u
    Y=n*T/X_l; % Inter-exceedances {T_i} normalized by the tail function (see, {Y_i} in  Algorithm 1, step 1)
    Z = sort(Y ); % Sorted normalised inter-exceedances
    L = length(Z);
    %L
   % Computation of a pilot intervals estimate (see, Algorithm 1, step 2) by formula (4) in Markovich and Rodionov (2022).   
                        max_theta=max(T);
    if (max_theta <=2 )
        theta(u_c)=min(1,(2*(sum(T))^2)/(length(find(T))*sum(T.^2)));
    else
        theta(u_c)=min(1,(2*(sum(T-1))^2)/((length(find(T))-1)*sum((T-1).*(T-2))));
    end
   k=floor(theta(u_c)*L); % The number of selected largest inter-cluster times (see, Algorithm 1, step 2)
  % k=floor(min(sqrt(L),theta(u_c)*L));
 % k=floor((log(L))^2);
   %k
   if (theta(u_c)==1)
       k=L-1;
   end
   if (L>1)
            Z_0=Z(L-k);
   else
       Z_0=0;
   end
  % Normalized W2 statistic (see, Algorithm 1, step 3).
   w_plus(u_c) = 1/(12*k);
     for j = 1 : k 
        
        w_plus(u_c) =w_plus(u_c)+ (1-exp(-theta(u_c)*(Z(j+L-k)-Z_0))-(j-0.5)/k).^2;
     end
                 if (L<40)
                w_plus(u_c) = (w_plus(u_c)-0.4/L+0.6/L.^2)*(1+1/L); 
                         else
                w_plus(u_c) = w_plus(u_c); 
             end
                                   
   %if(abs(w_plus(u_c) - sigma_1)<=eps) % Cheking of condition (13)
   if(w_plus(u_c)<= sigma_1)            % Cheking of condition (16)
   %if(w_plus(u_c)>= sigma_1)
       u_opt(u_c)=u_c;                  % Computation of threshold obtained as solution of discrepancy equation that satisfies (13) or (16).
       theta_opt(u_c)=theta(u_c);       % Computation of extremal index estimate with the selected threshold
          else
       f_m=f_m+1;
       u_opt(u_c)=0;
       theta_opt(u_c)=0;
   end
end
theta_opt
if(theta_opt==0)
    flag_s=1;
  else
       flag_s=0;
      b = theta_opt(theta_opt ~= 0);
   b
   c = u_opt(u_opt ~= 0);
   c
      c_min=min(c);
   c_min
   c_max=max(c);
   c_max
           l_m=length(u_opt)-f_m;
l_m
% Calculation of statistics by formula (15).
theta_1=mean(b);
theta_1
theta_2=theta_opt(c_min);
theta_2
theta_3=theta_opt(c_max);
theta_3
end
%CI=prctile(theta(:),[5,95]); % Confidence interval based on percintiles of  estimates
 
%[thr_val_q' theta']

%plot(thr_val_q, w_plus,thr_val_q,sigma_1)
%      xlabel('Quantile level')
% ylabel('w2')
     
%plot(thr_val_q, theta)
%     xlabel('Quantile level')
% ylabel('Extremal index')
 

end

