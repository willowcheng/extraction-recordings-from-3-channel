function output=MMSECohen2004(signal,fs,IS)

% Denoising using MMSE log-spectral Estimator
% output=MMSECohen2004(signal,fs,IS)
% signal is the noisy signal
% fs is the sampling frequency
% IS is the initial noise length in seconds (default 0.25 seconds)
% output is the denoised reconstructed signal
% functions required:
%       VAD
%       SEGMENT
%       OVERLAPADD
% This function uses the method proposed by I.Cohen 2004 for A Priori
% estimation: "SPEECH ENHANCEMENT USING A NONCAUSAL A PRIORI SNR ESTIMATOR"
% The MMSE gain is based on Ephraims (1985) method: "SPEECH ENHANCEMENT
% USING A MINIMUM MEAN SQUARE ERROR LOG-SPECTRAL AMPLITUDE ESTIMATOR"
%
% Esfandiar Zavarehei
% Dec 2004
% Last Modified: 24-01-05 (Ove


if (nargin<3 | isstruct(IS))
    IS=.25; %seconds
end
W=fix(.032*fs); %Window length is 25 ms

SP=.4; %Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
wnd=hamming(W);
if (nargin>=3 & isstruct(IS))%This option is for compatibility with another programme
    W=IS.windowsize
    SP=IS.shiftsize/W;
    
    wnd=IS.window;
    if isfield(IS,'IS')
        IS=IS.IS;
    else
        IS=.25;
    end
end

pre_emph=0;
signal=filter([1 -pre_emph],1,signal);
NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments

y=segment(signal,W,SP,wnd);
Y=fft(y);
YPhase=angle(Y(1:fix(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:fix(end/2)+1,:));%Specrogram
numberOfFrames=size(Y,2);
FreqResol=size(Y,1);
nfft=size(Y,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%         VAD Parameters
NoiseCounter=0;
NoiseLength=9;%This is a smoothing factor for the noise updating
%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=mean(Y(:,1:NIS)')'; %initial Noise Power Spectrum mean
LambdaD=mean((Y(:,1:NIS)').^2)';%initial Noise Power Spectrum variance

alpha=.98; %used in smoothing xi
G=ones(size(N));%Initial Gain used in calculation of the new xi
Gamma=G;
L=3;
b=[.25;.5;.25]; %The design is based on this length of b (e.g. calculation of Lph)
bGeneral=repmat(b,1,L+1); bGeneral(2,1)=0; % (n,i)~=(0,0) [eq. 14. Cohen 2004]
mu=.8;
xmin=10^(-25/10);
LambdaMin=xmin.*LambdaD; %Updated in VAD
Beta=2;

Gamma1p5=gamma(1.5); %Gamma function at 1.5


X=zeros(size(Y)); % Initialize X (make room)
h=waitbar(0,'Wait...');
for i=1:numberOfFrames
    %>>>>>>>>>>>>>>>>>> VAD AND NOISE UPDATE START
    if i<=NIS
        SpeechFlag=0;
        NoiseCounter=100;
    else
        [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i),N,NoiseCounter); %Magnitude Spectrum Distance VAD
    end
    if SpeechFlag==0
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1); %Update and smooth noise mean
        LambdaD=(NoiseLength*LambdaD+(Y(:,i).^2))./(1+NoiseLength); %Update and smooth noise variance
        LambdaMin=xmin.*LambdaD;
    end
    %<<<<<<<<<<<<<<<<<<  VAD AND NOISE UPDATE END

    %>>>>>>>>>>>>>>>>>>  Lphb Calculation [Eq. 14 Cohen 2004]     START
    LL=min(L,numberOfFrames-i);
    bbGeneral=bGeneral(:,1:LL+1);
    bGeneralNorm=sum(sum(bbGeneral));
    bSpecialNorm=sum(sum(bbGeneral(1:2,:)));
    for k=2:nfft-1 %Frequency constraints
        Lphb(k,1)=sum(sum(bbGeneral.*(Y(k-1:k+1,i:i+LL).^2)))/bGeneralNorm;
    end
    Lphb(1,1)=sum(sum(bbGeneral(2:3,:).*(Y(1:2,i:i+LL).^2)))/bSpecialNorm;
    Lphb(nfft,1)=sum(sum(bbGeneral(1:2,:).*(Y(nfft-1:nfft,i:i+LL).^2)))/bSpecialNorm;

    Lphb=Lphb-Beta*LambdaD;
    Lphb=max(Lphb,0);
    %<<<<<<<<<<<<<<<<<<  Lphb Calculation [Eq. 14 Cohen 2004]     END
    
    %>>>>>>>>>>>>>>>>>>  Lph Calculation [Eq. 13 Cohen 2004]      START
    if (i>1)
        for k=2:nfft-1
            Lph(k,1)=mu*(X(k,i-1).^2) + (1-mu)*(mu*sum(b.*Lh(k-1:k+1,1))...
                + (1-mu)*Lphb(k));
        end
        Lph(1,1)=mu*(X(k,i-1).^2) + (1-mu)*(mu*sum(b(2:3).*Lh(k:k+1,1))/sum(b(2:3))...
                + (1-mu)*Lphb(k));
        Lph(nfft,1)=mu*(X(k,i-1).^2) + (1-mu)*(mu*sum(b(1:2).*Lh(k-1:k,1))/sum(b(1:2))...
                + (1-mu)*Lphb(k));
    else %that is if i==1,
        % Initial conditions: |A(k,-1)|=0, AND, Lh(k,-1)=LambdaMin [Table
        % 1 Cohen 2004]
        Lph(:,1)=mu*LambdaMin + (1-mu)*Lphb;
    end
    %<<<<<<<<<<<<<<<<<<  Lph Calculation [Eq. 13 Cohen 2004]      END
    
    %>>>>>>>>>>>>>>>>>>  Lh Calculation [Eq. 12 Cohen 2004]      START
    Lh= (Lph./(Lph+LambdaD)) .* (LambdaD + ((Lph.*(Y(:,i).^2))./(Lph+LambdaD)));
    %<<<<<<<<<<<<<<<<<<  Lh Calculation [Eq. 12 Cohen 2004]      END
    
    xi=Lh./LambdaD;
    Gamma=(Y(:,i).^2)./LambdaD;
    nu=Gamma.*xi./(1+xi);

%     G=(Gamma1p5*sqrt(nu)./Gamma).*exp(-nu/2).*((1+nu).*besseli(0,nu/2) +nu.*besseli(1,nu/2)); %MMSE [Ephraim 1984] 
%     G(find(isnan(G)))=1;    
%     G(find(isinf(G)))=1;
    G= (xi./(1+xi)).*exp(.5*expint(nu)); % Log spectral MMSE [Ephraim 1985]
    X(:,i)=G.*Y(:,i);
    
    waitbar(i/numberOfFrames,h,num2str(fix(100*i/numberOfFrames)));
end
close(h);
output=OverlapAdd2(X,YPhase,W,SP*W);
output=filter(1,[1 -pre_emph],output);

function ReconstructedSignal=OverlapAdd2(XNEW,yphase,windowLen,ShiftLen);

%Y=OverlapAdd(X,A,W,S);
%Y is the signal reconstructed signal from its spectrogram. X is a matrix
%with each column being the fft of a segment of signal. A is the phase
%angle of the spectrum which should have the same dimension as X. if it is
%not given the phase angle of X is used which in the case of real values is
%zero (assuming that its the magnitude). W is the window length of time
%domain segments if not given the length is assumed to be twice as long as
%fft window length. S is the shift length of the segmentation process ( for
%example in the case of non overlapping signals it is equal to W and in the
%case of %50 overlap is equal to W/2. if not givven W/2 is used. Y is the
%reconstructed time domain signal.
%Sep-04
%Esfandiar Zavarehei

if nargin<2
    yphase=angle(XNEW);
end
if nargin<3
    windowLen=size(XNEW,1)*2;
end
if nargin<4
    ShiftLen=windowLen/2;
end
if fix(ShiftLen)~=ShiftLen
    ShiftLen=fix(ShiftLen);
    disp('The shift length have to be an integer as it is the number of samples.')
    disp(['shift length is fixed to ' num2str(ShiftLen)])
end

[FreqRes FrameNum]=size(XNEW);

Spec=XNEW.*exp(j*yphase);

if mod(windowLen,2) %if FreqResol is odd
    Spec=[Spec;flipud(conj(Spec(2:end,:)))];
else
    Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
end
sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);
weight=sig;
for i=1:FrameNum
    start=(i-1)*ShiftLen+1;    
    spec=Spec(:,i);
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+real(ifft(spec,windowLen));    
end
ReconstructedSignal=sig;

function Seg=segment(signal,W,SP,Window)

% SEGMENT chops a signal to overlapping windowed segments
% A= SEGMENT(X,W,SP,WIN) returns a matrix which its columns are segmented
% and windowed frames of the input one dimentional signal, X. W is the
% number of samples per window, default value W=256. SP is the shift
% percentage, default value SP=0.4. WIN is the window that is multiplied by
% each segment and its length should be W. the default window is hamming
% window.
% 06-Sep-04
% Esfandiar Zavarehei

if nargin<3
    SP=.4;
end
if nargin<2
    W=256;
end
if nargin<4
    Window=hamming(W);
end
Window=Window(:); %make it a column vector

L=length(signal);
SP=fix(W.*SP);
N=fix((L-W)/SP +1); %number of segments

Index=(repmat(1:W,N,1)+repmat((0:(N-1))'*SP,1,W))';
hw=repmat(Window,1,N);
Seg=signal(Index).*hw;

function [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(signal,noise,NoiseCounter,NoiseMargin,Hangover)

%[NOISEFLAG, SPEECHFLAG, NOISECOUNTER, DIST]=vad(SIGNAL,NOISE,NOISECOUNTER,NOISEMARGIN,HANGOVER)
%Spectral Distance Voice Activity Detector
%SIGNAL is the the current frames magnitude spectrum which is to labeld as
%noise or speech, NOISE is noise magnitude spectrum template (estimation),
%NOISECOUNTER is the number of imediate previous noise frames, NOISEMARGIN
%(default 3)is the spectral distance threshold. HANGOVER ( default 8 )is
%the number of noise segments after which the SPEECHFLAG is reset (goes to
%zero). NOISEFLAG is set to one if the the segment is labeld as noise
%NOISECOUNTER returns the number of previous noise segments, this value is
%reset (to zero) whenever a speech segment is detected. DIST is the
%spectral distance. 
%Saeed Vaseghi
%edited by Esfandiar Zavarehei
%Sep-04

if nargin<4
    NoiseMargin=3;
end
if nargin<5
    Hangover=8;
end
if nargin<3
    NoiseCounter=0;
end
    
FreqResol=length(signal);

SpectralDist= 20*(log10(signal)-log10(noise));
SpectralDist(find(SpectralDist<0))=0;

Dist=mean(SpectralDist); 
if (Dist < NoiseMargin) 
    NoiseFlag=1; 
    NoiseCounter=NoiseCounter+1;
else
    NoiseFlag=0;
    NoiseCounter=0;
end

% Detect noise only periods and attenuate the signal     
if (NoiseCounter > Hangover) 
    SpeechFlag=0;    
else 
    SpeechFlag=1; 
end 