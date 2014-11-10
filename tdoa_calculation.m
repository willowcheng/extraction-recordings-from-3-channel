%% Extraction a Target Speech Signal from 3-channel Recordings
% Copyright (c) 2014-2015, Ottawa-Carleton Institute for Electrical 
% and Computer Engineering of University of Ottawa
% Author: Liu Cheng & Hongyu Zou
% Student Number: 7486632, 
% Contact Email: lchen156@uottawa.ca

% This script is based on bss_locate_spec.m written by Charles Blandin 
% and Emmanuel Vincent in 2011, public example provided by MATLAB which
% is used as Phased Array Toolbox
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
d = 0.5;  
% read mixture audio
waves_mixture13 = [audioread('mixture1.wav'),audioread('mixture3.wav')];  
% calculate tdoa using parameters above
tdoa = bss_locate_spec(waves_mixture13, fs, d, nsrc,'GCC-PHAT');
% store values of 3 angles into vector theta
theta = [asind(tdoa(1)*c/d), asind(tdoa(2)*c/d), asind(tdoa(3)*c/d)];

%% Define a Uniform Linear Array
% in this step, using frost beamformer to extract directional source
% this is part of code provided in official example

% follow the range of frequency of human being hearing system
% to set microphone parameters
% ...
hmic = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 20e3]);

% ...
ha = phased.ULA(nmic,space,'Element',hmic);

% add elevation to meet requirement of system function
ang_speech = [theta; 0, 0, 0];

% choose second source
%...
angSteer = ang_speech(:,2);

nsamplePerframe = 8000;  % set buffer size for output audio

hsig = dsp.SignalSource('Signal',waves,...
    'SamplesPerFrame',nsamplePerframe);  % set source audio parameters

% set built-in audio player for playing audio
isAudioSupported = helperMicrophoneExampleAudioSupported;
if isAudioSupported   % make sure audio is supported 
    hap = dsp.AudioPlayer('SampleRate', fs);
end

% set parameters for beamformer, related information can be found in 
% official help documentation
hbf = phased.FrostBeamformer('SensorArray',ha,'SampleRate',fs,...
    'PropagationSpeed',c,'FilterLength',20,'DirectionSource','Input port');
reset(hsig);
out = zeros(nsample, 1);  % initial out as all zeros

% write generated samples frame by frame using given beamformer and audio
for m = 1:nsamplePerframe:nsample
    out(m:m+nsamplePerframe-1,:) = step(hbf,step(hsig),angSteer);
end

plot(t,out);  % plot time and waveform, set labels, title, limination of y
xlabel('Time (sec)'); ylabel ('Amplitude (V)');
title('Frost Beamformer Output'); ylim([-3 3]);
audiowrite('generate_frost2.wav', out, fs);  % write audio into audio file

%% Wiener Noise Filtering
% use wiener filter for noise deduction

wave = audioread('generate_frost2.wav');  % read audio file into wave
noisesamples = 2 * fs;  % number of noise samples

% function is used and modified as part of work of Liu Ming, 2008
[comparison_output, output] = WienerNoiseReduction(wave, fs, noisesamples);
% two methods are involved while Harmonic Regeneration Noise Reduction 
% method has better performance; thus choose it as output

audiowrite('output2.wav', output, fs);  % write audio into final audio wave