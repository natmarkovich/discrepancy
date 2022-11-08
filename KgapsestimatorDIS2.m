
function [theta_1, theta_2, theta_3, flag_s] = KgapsestimatorDIS2(X)
% K-gaps estimator by Sueveges and Davison (2010)
% Inter-cluster times T are replaced by K-gaps S with discrepancy method used to select u 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%K-gaps estimator (see, formula (6) in Markovich and Rodionov (2022))==============
X_l= length(X); % Length of the degree sequence we need for Interval Estimator
thr_val_q = [90 90.5 91 91.5 92 92.5 93 93.5 94 94.5 95 95.5 96 96.5 97 97.5 98 98.5 99 99.5]; % percentage quantiles
%thr_val_q = [95 95.5 96 96.5 97 97.5 98 98.5 99 99.5]; % percentage quantiles
%thr_val_q = [95 96  97  98  99];
q_1=20; % the number of quantiles
L=0; % the number of inter-cluster times
%s=0.05; % proportion of largest inter-cluster times
s=0.56;
eps=0.01; % upper deviation between w2 and the discrepancy value sigma_1
sigma_1=0.05; % The mode of the statistic w2 (see, formula (14) in Markovich and Rodionov (2022))
sigma_1=1.49; %  The 99.98% quantile of the statistic w2 used in formula (16) as sigma_2
theta_opt=0; % An initial optimal value of extremal index
f_m=0;       % Number of unseccesful solutions of the discrepancy equation
flag_s=0; % An indicator of a sample when there is no solution of the discrepancy equation 
theta_1=0; theta_2=0; theta_3=0;  %Initial values of statistics in formula (15) by Markovich and Rodionov (2022)
thetaK_opt=0;
for u_c=1:length(thr_val_q)
        u=prctile(X,thr_val_q(u_c)); % The treshold is equal to the quantile of level u_c
        %u
        flag=0; % An indicator of exceedances over threshold u
    %% Computation of inter-exceedance times T_i(u) (see, Algorithm 1, step 1)       
T_t=0; % Point process of exceedances
    T=0;Y=0; 
    n=0; % The number of inter-cluster times
  
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
    % n is the number of exceedances over threshold u
   % n
       U=0; KK=0; w_plus=0;
   for K=1:q_1
      S=max(T-K,0);
    Y=n*S/X_l; % K-gaps normalized by tail function
    Z = sort(Y);
               L = length(Y); % The number of K-gaps
   %% Computation of a pilot intervals estimate (see, Algorithm 1, step 2) by formula (4) in Markovich and Rodionov (2022).   
                   max_theta=max(T);
    if (max_theta <=2 )
        theta_1(u_c)=min(1,(2*(sum(T))^2)/(length(find(T))*sum(T.^2)));
    else
        theta_1(u_c)=min(1,(2*(sum(T-1))^2)/((length(find(T))-1)*sum((T-1).*(T-2))));
    end
  %  theta_1(u_c)
      %k=floor(min(sqrt(L),theta_1(u_c)*L));
      % k=floor((log(L))^2);
  k=floor(theta_1(u_c)*L);
   if (k>=L)
       k=L-1;
   end
          N_C=0; % The number of non-zero S_i (the number of clusters)
               for i=1:L
        if (Y(i) ~= 0)
        N_C=N_C+1;
        end
               end
        %% Computation of the K-gaps estimate (6) in Markovich and Rodionov (2022)
    N_CC(u_c,K)=N_C;
       a_1=n-1-N_C;
    a_2=2*N_C;
    a_3=sum(Y);
        theta(u_c,K)= 0.5*((a_1+a_2)/a_3+1-(((a_1+a_2)/a_3+1).^2-4*a_2/a_3).^0.5);
      theta(u_c,K)
          % Discrepancy method
     % Normalized W2 statistic (see, Algorithm 1, step 3).
      w_plus = 1/(12*N_C);
  %% Normalized W2 statistic (see, Algorithm 1, step 3).
          
           for j = 1 : k 
                    w_plus =w_plus+ (1-exp(-theta(u_c,K)*(Z(j+L-k)-Z(L-k)))-(j-0.5)/k).^2;
              end
    
                 if (L<40)
                w_plus = (w_plus-0.4/L+0.6/L.^2)*(1+1/L); 
             
             end
             
  %if(abs(w_plus - sigma_1)<=eps)  % Cheking of condition (13)
     if(w_plus <= sigma_1)         % Cheking of condition (16)
       u_opt(u_c,K)=u_c;           % Computation of pairs (threshold,K), solutions of discrepancy equation 
                                   % that satisfies (13) or (16) (see, Remark 4 in Markovich and Rodionov 2022).
       K_opt(u_c,K)=K;             
       theta_opt(u_c,K)=theta(u_c,K); % Computation of  K-gaps estimate of the extremal index with the selected threshold
      else
              f_m=f_m+1;
             u_opt(u_c,K)=0;
              K_opt(u_c,K)=0;
       theta_opt(u_c,K)=0;
      end
      %w_plus
   end
end
   


theta_opt

if(theta_opt==0)
    flag_s=1;
  else
       flag_s=0;
      b = theta_opt(theta_opt ~= 0);
   %b
   c = u_opt(u_opt ~= 0);
   %c
   d = K_opt(K_opt ~= 0);
   %d
      c_min=min(c);
   %c_min
   d_min=min(d);
   %d_min
   c_max=max(c);
   %c_max
   d_max=max(d);
   %d_max
        l_m=length(u_opt)-f_m;
%l_m
theta_1=mean(b);
theta_1

%theta_2=theta_opt(c_min,d_min);
theta_2a=theta_opt(c_min,:);
theta_2=mean(theta_2a);
theta_2a
%theta_3=theta_opt(c_max,d_max);
theta_3a=theta_opt(c_max,:);
theta_3=mean(theta_3a);
theta_3a
%plot(thr_val_q, theta, thr_val_q, theta_1,'b--o',thr_val_q, theta_2,'c*',thr_val_q, theta_3,'--')
plot(thr_val_q, s, thr_val_q, theta_1,'b--o',thr_val_q, theta_2,'c*',thr_val_q, theta_3,'--')
      xlabel('Quantile level')
 ylabel('Extremal index')
end

%if(thetaK_opt==0)
%    flag_s=1;
%  else
%       flag_s=0;
% end
%[thr_val_q' theta']

% plot(thr_val_q, w_plus)
%      xlabel('Quantile level')
% ylabel('w2')
     
 
 


end

