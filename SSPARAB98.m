function output=SSPARAB98(signal,fs,IS)

% output=SSPARAB98(signal,fs,IS)
% "A Parametric Formulation of the Generalized Spectral Subtraction Method"
% This is very similar to SSPARASIM98
% Spectral subtraction based on the above paper by Boh Lim Sim 1998. In
% their method an optimized estimator based on spectral subtraction
% assumptions is derived. This estimator uses estimates of a priori and a
% posteriori SNR in its gain function. Also a flooring function is
% incorporated which uses a function to floor very small values of
% amplitude estimate. the paper Hasan04 (Modified A Priori SNR for speech
% enhancemnet) presents this formulation as well which was alternatively
% used in implementation of this paper. Decision-Directed method is used
% for estimation of A priori SNR.
% Author: Esfandiar Zavarehei
% Created: MAR-05



if (nargin<3 | isstruct(IS))
    IS=.25; %Initial Silence or Noise Only part in seconds
end
W=fix(.025*fs); %Window length is 25 ms
SP=.4; %Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
wnd=hamming(W);

%IGNORE FROM HERE ...............................
if (nargin>=3 & isstruct(IS))%This option is for compatibility with another programme
    W=IS.windowsize
    SP=IS.shiftsize/W;
    %nfft=IS.nfft;
    wnd=IS.window;
    if isfield(IS,'IS')
        IS=IS.IS;
    else
        IS=.25;
    end
end
% ......................................UP TO HERE

pre_emph=0;
signal=filter([1 -pre_emph],1,signal);

NIS=fix((IS*fs-W)/(SP*W) +1);%number of initial silence segments

y=segment(signal,W,SP,wnd); % This function chops the signal into frames
Y=fft(y);
YPhase=angle(Y(1:fix(end/2)+1,:)); %Noisy Speech Phase
Y=abs(Y(1:fix(end/2)+1,:));%Specrogram
numberOfFrames=size(Y,2);
FreqResol=size(Y,1);

N=mean(Y(:,1:NIS)')'; %initial Noise Power Spectrum mean
LambdaD=mean((Y(:,1:NIS)').^2)';%initial Noise Power Spectrum variance
alpha=.99; %used in smoothing xi (For Deciesion Directed method for estimation of A Priori SNR)
PowExp=2; % Power Exponent
Beta=(gamma(PowExp+1)-(gamma(1+PowExp/2)^2))/gamma(PowExp+1);%[Sim98, Eq. 21]

NoiseCounter=0;
NoiseLength=9;%This is a smoothing factor for the noise updating
mu=.05;

X=zeros(size(Y)); % Initialize X (memory allocation)

h=waitbar(0,'Wait...');

for i=1:numberOfFrames
    %%%%%%%%%%%%%%%%VAD and Noise Estimation START
    if i<=NIS % If initial silence ignore VAD
        SpeechFlag=0;
        NoiseCounter=100;
    else % Else Do VAD
        [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i),N,NoiseCounter); %Magnitude Spectrum Distance VAD
    end
    
    if SpeechFlag==0 % If not Speech Update Noise Parameters
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1); %Update and smooth noise mean
        LambdaD=(NoiseLength*LambdaD+(Y(:,i).^2))./(1+NoiseLength); %Update and smooth noise variance
    end
    %%%%%%%%%%%%%%%%%%%VAD and Noise Estimation END
    
    Gamma=(Y(:,i).^2)./LambdaD; %A postiriori SNR
    if i>1
        xi=(alpha*(X(:,i-1).^2) + (1-alpha).*max(Y(:,i).^2 - LambdaD,0))./LambdaD;
    else
        xi=(ones(size(Gamma))+Gamma)/2;
    end

    X(:,i)=(((xi.^PowExp)./((xi.^PowExp)+Beta)).*(Y(:,i).^PowExp - N.^PowExp)).^(1/PowExp);
    if i>1
        X(:,i)=max(X(:,i),.5.*(mu.*Y(:,i)+X(:,i-1)));
    end
        
            
    waitbar(i/numberOfFrames,h,num2str(fix(100*i/numberOfFrames)));
end

close(h);
output=OverlapAdd2(X,YPhase,W,SP*W); %Overlap-add Synthesis of speech
output=filter(1,[1 -pre_emph],output); %Undo the effect of Pre-emphasis

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