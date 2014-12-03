%% Extraction a Target Speech Signal from 3-channel Recordings
% Copyright (c) 2014-2015, Ottawa-Carleton Institute for Electrical 
% and Computer Engineering of University of Ottawa
% Author: Liu Cheng & Hongyu Zou
% Student Number: 7486632 & 7493642 
% Contact Email: lchen156@uottawa.ca & hzou047@uottawa.ca

% This script is based on bss_locate_spec.m written by Charles Blandin 
% and Emmanuel Vincent in 2011
% All work was modified by Liu Cheng & Hongyu Zou in Nov 2014 to 
% apply related algorithms for actual application 

%% General parameters setting
clc;
clear all;
close all;

fs = 16000;  % sample frequency is set to 16kHz
c = 343;  % propagation rate of sound in air (m/s)
nsrc = 3;  % indicate number of source of three
nmic = 3;  % similar, indicate number of microphones
nsample = fs * 12;  % total samples in all three mixture audio
space = 0.25;  % space for each micphones is 25cm
t = 0:1/fs:(12-1/fs); % total time length in terms of samples
% read three mixture audio into waves
waves = [audioread('mixture1.wav'), audioread('mixture2.wav'),...
    audioread('mixture3.wav')];


%% Calculate direction of 3 audio sources and generalize parameters
% Thanks for the work of Charles Blandin and Emmanuel Vincent, related
% MATLAB scripts are used and modified for calculation
% 'stft_multi.m' is used as dependent function for bss_locate_spec function
% by work of Charles Blandin and Emmanuel Vincent

% distance between microphone 1 and microphone 3
% for applying calculation of TDOA
d = 0.25;  
% read mixture audio
waves_mixture23 = [audioread('mixture2.wav'),audioread('mixture3.wav')];  
% calculate tdoa using parameters above
tdoa = bss_locate_spec(waves_mixture23, fs, d, nsrc,'GCC-PHAT');
% store values of 3 angles into vector theta
theta = [asind(tdoa(1)*c/d), asind(tdoa(2)*c/d), asind(tdoa(3)*c/d)];

%% Beamforming Filtering
% Start by scratch£¬ followed by equation of realization
M=3; % Size of Rx Matrix
space=0.25; % distance between two sensors
t0=space/c;
a=0.985; % a parameter that controls the level of the point source noise relative to that of the spatially white noise
sita0=90-theta(3); % The direction of desired signal
sita1=90-theta(1); % The direction of interference1
sita2=90-theta(2); % The direction of interference2
y1=wavread('mixture1.wav');
y2=wavread('mixture2.wav');
y3=wavread('mixture3.wav');
Y1=fftshift(fft(y1'));
Y2=fftshift(fft(y2'));
Y3=fftshift(fft(y3'));
V=[Y1;Y2;Y3];
k=1;
I=eye(3); % white noise
signal=sita0*pi/180; 
IR1=sita1*pi/180;    
IR2=sita2*pi/180;    
for f=-fs/2:fs/(192000-1):fs/2
    v=[1;exp(-i*t0*2*pi*f*cos(signal));exp(-i*2*t0*2*pi*f*cos(signal))]; % steer vector of desired signal
    I1=[1;exp(-i*t0*2*pi*f*cos(IR1));exp(-i*2*t0*2*pi*f*cos(IR1))]; % steer vector of interference1
    I2=[1;exp(-i*t0*2*pi*f*cos(IR2));exp(-i*2*t0*2*pi*f*cos(IR2))]; % steer vector of interference2
    Sn=(1-a)*I+a*(I1*I1'+I2*I2'); % psn
    WW=v'*inv(Sn)*v;                         
    W=inv(Sn)*v/WW;
    H=W';
    B(k)=H*V(:,k);
    k=k+1;  
end
Z=[B(96001:192000),B(1:192000)];
X=ifft(Z,192000);
x=real(X);
wavwrite(x,16000,'Original_0.985.wav');

%% Parametric Spectral Subtraction
% directly use SSPARAB98 function which implement Parametric Spectral
% Subtraction method
output=SSPARAB98(x,fs,2);
wavwrite(output,16000,'output.wav');

% output1=SSBoll79(x,fs,2);
% wavwrite(output1,16000,'output1.wav');
% 
% output2=MMSECohen2004(x,fs,2);
% wavwrite(output2,16000,'output2.wav');
% 
% output3=MMSESTSA85(x,fs,2);
% wavwrite(output3,16000,'output3.wav');
% 
% output4=MMSESTSA84(x,fs,2);
% wavwrite(output4,16000,'output4.wav');
% 
% output5=WienerScalart96(x,fs,2);
% wavwrite(output5,16000,'output5.wav');
% 
% output6=SSBerouti79(x,fs,2);
% wavwrite(output6,16000,'output6.wav');
% 
% output7=SSMultibandKamath02(x,fs,2);
% wavwrite(output7,16000,'output7.wav');
% 
% output8=SSScalart96(x,fs,2);
% wavwrite(output8,16000,'output8.wav');
% 
% output9=SSPARAB98(x,fs,2);
% wavwrite(output9,16000,'output9.wav');
