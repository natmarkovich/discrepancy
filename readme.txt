Description of files


1. Main file main_1 

               implements Algorithm 1 by Markovich & Rodionov (2022) Threshold selection for extremal index estimation.
arXiv:2009.02318v1
to simulated data where the intervals estimator is used to estimate the extremal index. By Algorithm 1  threshold u is found among 
high quantiles of an underlying sample X that can be solutions of the discrepancy equation. 

Function generation_X(X(:,j)) 

               generates samples from several stochastic processes including the processes described in Section 4.1 by M & R (2022).

X is a matrix that includes all generated samples.

j=1,2,...,sample is a number of the sample, where sample denotes the number of repetitions.

The function generation_X returns back values of the following variables (theta_1, theta_2, theta_3, flag_s)

Function intervalsestimator(X) 

               implements the computation of the intervals estimator (see, (4) in M & R (2022)) and checking 
the fulfillment of inequalities (13) or (16) in M & R (2022) for high quantiles. The  number k of the largest order statistics 
may be taken as indicated in Algorithm 1, Step (2).

X is the jth sample generated.

The function intervalsestimator returns back values of the following variables (theta_1, theta_2, theta_3, flag_s). Here, theta_1, theta_2, theat_3 are values calculated by formula (15)  by M & R (2022).
flag_s indicates whether 
a solution of the discrepancy equation satisfying to inequality (13) or (16)  by M & R (2022) is found for the underlying sample. If flag_s=0, then the solution has been found, otherwise flag_s=1. 

2. Main file mainA1 

Function intervalsestimatorA1(X)

               implement Algorithm 1 by M & R (2022) where the threshold u for the intervals estimate is selected 
by the "plateau-finding" algorithm A1 by 
Ferreira (2018a) Heuristic tools for the estimation of the extremal index: a comparison of methods.
REVSTAT -- Statistical Journal, 16(1), 115--136.

X is the jth sample generated.

The function intervalsestimatorA1 returns back values of the following variables (theta_1,  flag_s). Here, theta_1 denotes the estimate obtained by algorithm A1 by 
Ferreira (2018a).
 flag_s indicates whether the latter estimate is found or not.

3. Main file ReadDataDewPointTemperatures 

               implements Algorithm 1 by M & R (2022) to real data set dewpoint_data2, see Sec. 5.2 and Tab.4 
by M & R (2022) for the description and results.

Function KgapsestimatorDIS2(X) 

                implements the K-gaps estimator (see (6) in M & R (2022)) where the threshold u is selected 
by the discrepancy method. The latter estimator is denoted as Kdis in Tab.4 by M & R (2022).

X is the jth sample generated. 

The function KgapsestimatorDIS2 returns back values of the following variables (theta_1, theta_2, theta_3, flag_s). Here, theta_1, theta_2, theat_3 and flag_s mean the same as in Item 1.
