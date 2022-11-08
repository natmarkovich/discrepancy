function X = generation_X(distribution_type, n, sample, MM_parametr,...
AR_parametr, ARMAX_parametr, AR2_parametr,MA2_parametr,GARCH_parametr,LINDLEY_parametr)
%% generation_X generates  matrix(n,sample) of a sample of random variables of simulated processes
%  lognormal 1,weibull 2, frechett 3, AR 4
% The desciption of the processes is given in Sec. 4.1 in Markovich and
% Rodionov (2022)


  %MM process
    if(distribution_type == 1)
        gener_frechett = @(t)(-log(t)).^(- 1);
                       for j = 1 : sample
                           nn=n+MM_parametr.m;
                           Z = gener_frechett(rand(nn, 1));
                                                     
            %Z = gener_frechett(rand(n+MM_parametr.m, 1));
                      % T = number of observations
            T = n+MM_parametr.m;
            t_1 = 1+MM_parametr.m;
                                    for t = t_1 : T
                m_1= max(MM_parametr.a1.*Z(t), MM_parametr.a2.*Z(t-1));
                %m_1
                m_2= max(MM_parametr.a3.*Z(t-2),MM_parametr.a4.*Z(t-3));
                %m_2
                X(t-MM_parametr.m, j) = max(m_1,m_2);
                %MM_parametr.a3.*Z(t-2,j),MM_parametr.a4.*Z(t-3,j));
            end
        end
    end
             
      % AR process with uniform noise
    if(distribution_type == 2)    
        for j = 1 : sample
            % T = number of observations
            T = n;
            % r = coefficient 
            X(1, j) =  rand(1, 1);
                        f = floor((AR_parametr.r - 1)* rand(n, 1)+0.5); 
            
                                   epsilon = f/AR_parametr.r;
                                    %epsilon = (AR_parametr.r - 1)/AR_parametr.r.* rand(n, 1);
             %epsilon
                        for t = 2 : T
                X(t, j) = 1/AR_parametr.r * X(t - 1, j) + epsilon(t, 1);
            end
        end
    end
    
    %ARMAX process
    if(distribution_type == 3)
        gener_frechett = @(t)(-log(t)).^(- 1);        
        for j = 1 : sample
            Z = gener_frechett(rand(n, 1));
            % T = number of observations
            T = n;
            % r = coefficient 
            X(1, j) =  Z(1,1);
            for t = 2 : T
                X(t, j) = max(ARMAX_parametr.alfa.*X(t-1, j),...
                    (1 - ARMAX_parametr.alfa).*Z(t));
            end
        end
    end
    
    % AR process with Cauchy noise
    if(distribution_type == 4)    
        for j = 1 : sample
            % T = number of observations
            T = n;
            % r = coefficient 
            X(1, j) =  rand(1, 1);
            %epsilon = stat::cauchyRandom(0,1,n);
            %rng('default');  % For reproducibility
            %X(1, j) =  trnd(1,1, 1);
            rng('default');  % For reproducibility
epsilon = trnd(1,n,1);
                        for t = 2 : T
                X(t, j) =0.7 * X(t - 1, j) + epsilon(t, 1);
            end
        end
    end
    %
     % AR(2) process with Pareto noise
    if(distribution_type == 5) 
        gener_pareto = @(t)(1-t.^(- 1/AR2_parametr.alpha));
        for j = 1 : sample
            epsilon = gener_pareto(rand(n, 1));
            % T = number of observations
            T = n;
            % r = coefficient 
            % rng('default');
            %X(1, j) =  gprnd(2,1,0);
            % rng('default');
            %X(2, j) =  gprnd(2,1,0);
            X(1, j) = 0; X(2, j) = 0;
                       % rng('default');  % For reproducibility
%epsilon = gprnd(2,1,0,n,1); % Pareto with tail index =2 
                        for t = 3 : T
                X(t, j) =0.95 * X(t - 1, j)-0.89*X(t - 2, j) + epsilon(t, 1);
            end
        end
    end
     % MA(2) process with Pareto noise
    if(distribution_type == 6) 
        gener_pareto = @(t)(t.^(- 1/MA2_parametr.alpha));
        for j = 1 : sample
            Z = gener_pareto(rand(n+2, 1));
            % T = number of observations
            T = n;            
                        for t = 3 : T+2
                X(t-2, j) =MA2_parametr.p * Z(t - 2)+MA2_parametr.q*Z(t - 1) + Z(t);
            end
        end
    end
    % GARCH(1,1) process with Gaussian noise
    if(distribution_type == 7) 
                for j = 1 : sample
                       % T = number of observations
            T = n; 
            Z = normrnd(0,1,[n 1]); 
            X(1, j) = 0; s(1, j) = 0;
                        for t = 2 : T
s(t, j) =(GARCH_parametr.alpha+GARCH_parametr.lambda * X(t-1,j).^2+GARCH_parametr.beta*s(t - 1,j).^2).^0.5;
% s(t, j)
% Z(t)
X(t, j) = s(t, j)*Z(t);             
                            end
        end
    end
    % Set by Marta Ferreira (2018)
    %MM with m=2
    if(distribution_type == 8)
        gener_frechett = @(t)(-log(t)).^(- 1);
                       for j = 1 : sample
                           nn=n+MMm2_parametr.m;
                           Z = gener_frechett(rand(nn, 1));
                                                     
            %Z = gener_frechett(rand(n+MM_parametr.m, 1));
                      % T = number of observations
            T = n+MMm2_parametr.m;
            t_1 = 1+MMm2_parametr.m;
                                    for t = t_1 : T
                m_1= max(MMm2_parametr.a1.*Z(t), MMm2_parametr.a2.*Z(t-1));
                m_2= MMm2_parametr.a3.*Z(t-2);
                X(t-MMm2_parametr.m, j) = max(m_1,m_2);
                %MM_parametr.a3.*Z(t-2,j),MM_parametr.a4.*Z(t-3,j));
            end
        end
    end
    
     % AR process with uniform noise with negative correlation
    if(distribution_type == 9)    
        for j = 1 : sample
            % T = number of observations
            T = n;
            % r = coefficient 
             X(1, j) =  rand(1, 1);
            % if(AR_parametr.r==2)
            f = ceil((AR_parametr.r - 1)* rand(n, 1)+0.5); 
           % f
           % else
           % f = ceil((AR_parametr.r - 1)* rand(n, 1));
           % f
           % end
             epsilon = f/AR_parametr.r;
                        for t = 2 : T
                X(t, j) = -1/AR_parametr.r * X(t - 1, j) + epsilon(t, 1);
            end
        end
    end
     %Lindley process
    if(distribution_type == 10)
        % Lindley_parametr.lambda, Lindley_parametr.mu - arrival and
        % service rates is doble-exponential or Laplace distributed
        %gener_doubleexponential = @(t)(-log((1+LINDLEY_parametr.mu/LINDLEY_parametr.lambda)*(1-t))./LINDLEY_parametr.mu); 
       % T = number of observations
            T = n;
         
            w = (LINDLEY_parametr.lambda+LINDLEY_parametr.mu)./LINDLEY_parametr.lambda;
            y = 1/LINDLEY_parametr.mu;
            
                        u = rand(T,1);
                        
                    Z = 0;
                    Z(u<=0.5) = y*log(w*u(u<=0.5));
                    Z(u>0.5) =  - y*log(w*(1-u(u>0.5)));
             
        for j = 1 : sample
         %   Z = gener_doublyexponential(rand(n, 1));

                        X(1, j) =  0;
            for t = 2 : T
                
              %  Z(t-1)
              %  X(t-1, j)
                X(t, j) = max(0,X(t-1, j)+Z(t-1));
               % X(t, j)
            end
           %  plot(X)
        end
        
       
 %     xlabel('t')
 %ylabel('Lindley process, X')
    end
    
   % u = rand( sampleSize );
   %                 out = zeros( sampleSize );
   %                 out(u<=0.5) = theta + lam*log(2*u(u<=0.5));
   %                 out(u>0.5) = theta - lam*log(2*(1-u(u>0.5)));
end

