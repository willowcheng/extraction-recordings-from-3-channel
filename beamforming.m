clear all; 
close all;
clc
M=3; % Size of Rx Matrix
d=0.25;
c=343;
t0=d/c;
fs=16000;
a=0.985;
sita0=90-10.2403; %信号方向
sita1=90-81.4509; %干扰方向1
sita2=90-(-81.4509); %干扰方向2
y1=wavread('mixture1.wav');
y2=wavread('mixture2.wav');
y3=wavread('mixture3.wav');
Y1=fftshift(fft(y1'));
Y2=fftshift(fft(y2'));
Y3=fftshift(fft(y3'));
V=[Y1;Y2;Y3];
k=1;
I=eye(3);

signal=sita0*pi/180; %信号的入射角度
IR1=sita1*pi/180;    %干扰信号入射方向
IR2=sita2*pi/180;    %干扰信号入射方向

for f=-fs/2:fs/(192000-1):fs/2
    v=[1;exp(-i*t0*2*pi*f*cos(signal));exp(-i*2*t0*2*pi*f*cos(signal))];
    I1=[1;exp(-i*t0*2*pi*f*cos(IR1));exp(-i*2*t0*2*pi*f*cos(IR1))];
    I2=[1;exp(-i*t0*2*pi*f*cos(IR2));exp(-i*2*t0*2*pi*f*cos(IR2))];
    Sn=(1-a)*I+a*(I1*I1'+I2*I2'); 
    WW=v'*inv(Sn)*v;                         
    W=inv(Sn)*v/WW;
    H=W';
    B(k)=H*V(:,k);
    k=k+1;  
end
Z=[B(96001:192000),B(1:192000)];
X=ifft(Z,192000);

% X1=ifftshift(ifft(B,192000));
% X2=ifft(B,192000);
x=real(X);
wavwrite(x,16000,'new_0.985.wav');


% theta=-90:1:90;%ULA估计角度变化的范围和频率选择
% V(:,:)=exp(-j*pi*K'*sin(theta*pi/180));   %20*1
%     %Sn=V*V'; %没有干扰时为此式 噪声谱矩
% B=W'*V;  

%1*1
% B(isnan(B))=1;
% x=ifft(B);
% x=real(X);
% sound(X,16000);


% P=abs(B);
% P = 20*log10(P/max(P));
% plot(theta,P);
% grid on;
