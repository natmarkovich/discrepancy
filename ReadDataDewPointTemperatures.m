clc;
clear;
%% Extremal index estimates for dewpoint temperatures data (see, Sec 5.2 in Markovich and Rodionov 2022)
%% Dewpoint temperatures data by Raymond et al. (2020)
M=dlmread('dewpoint_data2.txt');
s_1=M(:,1); s_2=M(:,2);s_3=M(:,3); s_4=M(:,4); s_5=M(:,5);s_6=M(:,6); 
s_7=M(:,7);s_8=M(:,8); s_9=M(:,9); s_10=M(:,10); s_11=M(:,11); s_12=M(:,12);
Z_l= length(s_1); % Length of the sample we use to estimate extremal index
Z_l
z_1=s_1(isfinite(s_1(:)), :);
min(z_1)
max(z_1)
length(z_1)
for j=1:length(z_1)
    thr_val_q(j)=j;
end
plot(thr_val_q, z_1,'k')
      xlabel('i')
 ylabel('daily-maximum dewpoint temperatures')
[theta_1, theta_2, theta_3, flag_s] = KgapsestimatorDIS2(z_1);
flag_s
theta_1
theta_2
theta_3