function [theta_1, flag_s] = intervalsestimatorA1(X)
% Intervals estimator with plateau-finding algorithm A1 Marta Ferreira
% (2018)
%   Detailed explanation goes here
%IE estimator==============================================
X_l= length(X); % Length of the degree sequence we need for Interval Estimator
thr_val_q = [90 90.5 91 91.5 92 92.5 93 93.5 94 94.5 95 95.5 96 96.5 97 97.5 98 98.5 99 99.5]; % percentage quantiles
L=0; % the number of inter-cluster times
%w=0.005;
%w=0.3;
w=0.25;
l=length(thr_val_q); % The number of intervals estimates
l
d_1=ceil(w*l); % The window size of moving average
d_1
m=ceil((l-2*d_1).^0.5); % PLateau length
m
theta_1=0; 
for u_c=1:length(thr_val_q)
        u=prctile(X,thr_val_q(u_c)); % The threshold selection as quantiles of levels thr_val_q
%% Computation of inter-exceedance times T_i(u) (see, Algorithm 1, step 1)       

        flag=0;
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
    % n is the number of exceedances over threshold u
    Y=n*T/X_l; % Inter-exceedances {T_i} normalized by the tail function (see, {Y_i} in  Algorithm 1, step 1)
    Z = sort(Y );  % Sorted normalised inter-exceedances
    L = length(Z); % The number of inter-exceedance times
    L
    % Computation of a pilot intervals estimate (see, Algorithm 1, step 2) by formula (4) in Markovich and Rodionov
                            max_theta=max(T);
    if (max_theta <=2 )
        theta(u_c)=min(1,(2*(sum(T))^2)/(length(find(T))*sum(T.^2))); % Intervals estimator
    else
        theta(u_c)=min(1,(2*(sum(T-1))^2)/((length(find(T))-1)*sum((T-1).*(T-2))));
    end
 
end
theta
s=std(theta); % The empirical standard deviation of intervals estimates for each threshold u
s

for k=1:l-2*d_1-1
    b_1=theta(k:k+2*d_1+1);
    b_1
    b(k)=mean(b_1);% The mean over a moving window of the size d_1
    b(k)
    k
end
c=0;
for j_1=1:l-2*d_1-m-1 % j_1 is the number of plateau
    j_1
    c(j_1)=0;
    for i_1=j_1+1:j_1+m-1
        i_1
    c(j_1)=c(j_1)+abs(b(i_1)-b(j_1));% The sum of deviations within one plateau
    end
    c(j_1)
    s
    2*s
    if(c(j_1)<=2*s)
        b_2=b(j_1:j_1+m-1);
     theta_opt(j_1)=mean(b_2); 
    else
       theta_opt(j_1)=0; 
    end
end
theta_opt
if(theta_opt==0)
    flag_s=1;
  else
       flag_s=0;
      b = theta_opt(theta_opt ~= 0);
   b
  
theta_1=b(1);
theta_1


end
%[thr_val_q' theta']

%plot(thr_val_q, w_plus,thr_val_q,sigma_1)
%      xlabel('Quantile level')
% ylabel('w2')
     
%plot(thr_val_q, theta)
%     xlabel('Quantile level')
% ylabel('Extremal index')
 
 

end
