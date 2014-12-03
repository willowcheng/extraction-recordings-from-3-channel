clear all; 
close all;
clc
M=3; % Size of Rx Matrix 
d=0.25;


K=0:M-1;
signal=10*pi/180; %信号的入射角度
v=exp(-j*K'*pi*sin(signal)) ; %阵列流形,信号源决定行数，阵元数决定列数   20*1
IR=-10/180*pi;    %干扰信号入射方向
I=exp(-j*K'*pi*sin(IR)) ; %干扰信号阵列流形

theta=-90:1:90;%ULA估计角度变化的范围和频率选择

    V(:,:)=exp(-j*pi*K'*sin(theta*pi/180));   %20*1
    %Sn=V*V'; %没有干扰时为此式 噪声谱矩阵
    Sn=eye(M)+I*I';  % 加入干扰后的噪声谱矩阵，《最优阵列处理》333页，式（6.76）
    WW=v'*inv(Sn)*v;                         %1*1
    W=inv(Sn)*v/WW;                    %20*1
    B=W'*V;                %1*1
    P=abs(B);

P = 20*log10(P/max(P));
plot(theta,P);
grid on;
