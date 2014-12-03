    function [tau, spec] = bss_locate_spec(x, fs, d, nsrc, local, pooling, tau_grid)

% BSS_LOCATE_SPEC Estimation of the source TDOAs in a stereo convolutive
% mixture using an angular spectrum
%
% [tau, spec] = bss_locate_spec(x, fs, d, nsrc, local, pooling, tau_grid)
%
% Inputs:
% x: nsampl x 2 matrix containing a stereo mixture signal
% fs: sampling frequency in Hz
% d: microphone spacing in meters
% nsrc: number of sources
% local: local angular spectrum function: 'GCC-PHAT', 'GCC-NONLIN'
%     (default), 'MUSIC', 'DS', 'MVDR', 'DNM', 'DSW' or 'MVDRW'
% pooling: pooling function: 'max' (default) or 'sum'
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds (default: 181
%     values linearly spaced between -d/343 and d/343)
%
% Outputs:
% tau: 1 x nsrc vector of estimated TDOAs in seconds
% spec: ngrid x 1 angular spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2010-2011 Charles Blandin and Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Charles Blandin, Emmanuel Vincent and Alexey Ozerov, "Multi-source TDOA
% estimation in reverberant audio using angular spectra and clustering",
% Signal Processing 92, pp. 1950-1960, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Errors and default values %%%
if nargin<4, error('Not enough input arguments.'); end
[nsampl,nchan]=size(x);
if nchan>nsampl, error('The input signal must be in columns.'); end
if nchan~=2, error('The input signal must be stereo.'); end
if nargin < 5, local = 'GCC-NONLIN'; end
if ~any(strcmp(local, {'GCC-PHAT' 'GCC-NONLIN' 'MUSIC' 'DS' 'MVDR' 'DNM' 'DSW' 'MVDRW'})), error('Unknown local angular spectrum.'); end
if nargin < 6, pooling = 'max'; end
if ~any(strcmp(pooling, {'max' 'sum'})), error('Unknown pooling function.'); end
if nargin < 7, tau_grid = linspace(-d/343, d/343, 181); end

%%% Time-frequency transform %%%
wlen = 1024;
if strfind(local, 'GCC'),
    % Linear transform
    X = stft_multi(x.',wlen);
    X = X(2:end,:,:);
else   
    % Quadratic transform
    hatRxx = qstft(x,fs,wlen,8,2);
    hatRxx = permute(hatRxx(:,:,2:end,:),[3 4 1 2]);
end
f = fs/wlen*(1:wlen/2).';

%%% Computing the angular spectrum %%%
% Computing the local angular spectrum
switch local
    case 'GCC-PHAT'
        spec = phat_spec(X, f, tau_grid);
    case 'GCC-NONLIN'
        c = 343;
        alpha = 10*c/(d*fs);
        spec = nonlin_spec(X, f, alpha, tau_grid);
    case 'MUSIC'
        spec = music_spec(hatRxx, f, tau_grid);
    case 'DS'
        spec = ds_spec(hatRxx, f, tau_grid);
    case 'MVDR'
        spec = mvdr_spec(hatRxx, f, tau_grid);
    case 'DNM'
        spec = dnm_spec(hatRxx, f, d, tau_grid);
    case 'DSW'
        spec = dsw_spec(hatRxx, f, d, tau_grid);
    case 'MVDRW'
        spec = mvdrw_spec(hatRxx, f, d, tau_grid);
end
% Pooling
switch pooling
    case 'max'
        spec = shiftdim(max(sum(spec,1),[],2));
    case 'sum'
        spec = shiftdim(sum(sum(spec,1),2));
end

%%% Finding the peaks %%%
[peaks, ind] = findpeaks(spec, 'minpeakdistance',2, 'sortstr','descend');
tau = tau_grid(ind(1:min(length(ind),nsrc)));

return;



function spec = phat_spec(X, f, tau_grid)

% PHAT_SPEC Computes the GCC-PHAT spectrum as defined in
% C. Knapp, G. Carter, "The generalized cross-correlation method for
% estimation of time delay", IEEE Transactions on Acoustics, Speech and
% Signal Processing, 24(4):320â€?27, 1976.
%
% spec = phat_spec(X, f, tau_grid)
%
% Inputs:
% X: nbin x nfram x 2 matrix containing the STFT coefficients of the input
%     signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of angular spectrum values

[nbin,nfram] = size(X(:,:,1));
ngrid = length(tau_grid);
X1 = X(:,:,1);
X2 = X(:,:,2);

spec = zeros(nbin,nfram,ngrid);
P = X1.*conj(X2);
P = P./abs(P);
for ind = 1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    spec(:,:,ind) = real(P.*EXP);
end

return;



function spec = nonlin_spec(X, f, alpha, tau_grid)

% NONLIN_SPEC Computes the nonlinear GCC-PHAT spectrum defined in
% B. Loesch, B. Yang, "Blind source separation based on time-frequency
% sparseness in the presence of spatial aliasing", in 9th Int. Conf. on
% Latent Variable Analysis and Signal Separation (LVA/ICA), pp. 1â€?, 2010.
%
% spec = nonlin_spec(X, f, alpha, tau_grid)
%
% Inputs:
% X: nbin x nfram x 2 matrix containing the STFT coefficients of the input
%     signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% alpha: nonlinearity parameter
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of angular spectrum values

[nbin,nfram] = size(X(:,:,1));
ngrid = length(tau_grid);
X1 = X(:,:,1);
X2 = X(:,:,2);

spec = zeros(nbin,nfram,ngrid);
P = X1.*conj(X2);
P = P./abs(P);
for ind = 1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    spec(:,:,ind) = ones(nbin,nfram) - tanh(alpha*real(sqrt(2-2*real(P.*EXP))));
end

return;



function spec = music_spec(hatRxx, f, tau_grid)

% MUSIC_SPEC Computes the MUSIC spectrum as defined in
% R. Schmidt, "Multiple emitter location and signal parameter estimation",
% IEEE Transactions on Antennas and Propagation, 34(3):276â€?80, 1986.
%
% spec = music_spec(hatRxx, f, tau_grid)
%
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values

[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = real(hatRxx(:,:,1,1));
R12 = hatRxx(:,:,1,2);
R22 = real(hatRxx(:,:,2,2));
TR = R11 + R22;
DET = R11.*R22 - abs(R12).^2;
lambda = .5*(TR + sqrt(TR.^2 - 4*DET));
V2 = (lambda-R11)./R12;

spec = zeros(nbin,nfram,ngrid);
for ind = 1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    spec(:,:,ind) = 1 ./ (1 - .5 * abs(1 + V2 .* conj(EXP)).^2./(1 + abs(V2).^2));
end

return;



function spec = ds_spec(hatRxx, f, tau_grid)

% DS_SPEC Computes the SNR in all directions using the DS beamformer
%
% spec = ds_spec(hatRxx, f, tau_grid)
%
% Inputs:
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values

[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);

SNR = zeros(nbin,nfram,ngrid);
for ind=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    SNR(:,:,ind) = (TR + 2*real(R12.*EXP))./(TR - 2*real(R12.*EXP));
end
spec = SNR;

return;



function spec = mvdr_spec(hatRxx, f, tau_grid)

% MVDR_SPEC Computes the SNR in all directions using the MVDR beamformer
%
% spec = mvdr_spec(hatRxx, f, tau_grid)
%
% Inputs:
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values

[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);

SNR = zeros(nbin,nfram,ngrid);
for ind=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    NUM = real(R11.*R22 - R12.*R21)./(TR - 2*real(R12.*EXP));
    SNR(:,:,ind) = NUM./(.5*TR-NUM);
end
spec = SNR;
        
return;



function spec = dnm_spec(hatRxx, f, d, tau_grid)

% APR_SPEC Computes the SNR in all directions using ML under a diffuse
% noise model
%
% spec = dnm_spec(hatRxx, f, d, tau_grid)
%
% Inputs:
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% d: microphone spacing in meters
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values

[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = real(hatRxx(:,:,1,1));
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = real(hatRxx(:,:,2,2));
c = 343;
SINC = sinc(2*f*d/c);
SINC2 = SINC.^2;

% Initializing the variances
vs = zeros(nbin,nfram,ngrid);
vb = zeros(nbin,nfram,ngrid);
for ind = 1:ngrid,
    
    % Computing inv(A) = [invA11 invA12; conj(invA11) -invA12]
    EXP = exp(-2*1i*pi*tau_grid(ind)*f);
    P = SINC .* EXP;
    invA11 = sqrt(.5)./(1-real(P)).*(1-conj(P));
    invA12 = -(1-P)./(SINC-EXP).*invA11;
    
    % Computing inv(Lambda) = [.5 invL12; 0 invL22]
    DEN = .5./(1-2*real(P)+SINC2);
    invL12 = (SINC2-1).*DEN;
    invL22 = 2*(1-real(P)).*DEN;
    
    % Computing vs and vb without nonnegativity constraint
    ARA1 = repmat(abs(invA11).^2,1,nfram).*R11 + repmat(abs(invA12).^2,1,nfram).*R22;
    ARA2 = ARA1 - 2 * real(repmat(invA11.*invA12,1,nfram).*R21);
    ARA1 = ARA1 + 2 * real(repmat(invA11.*conj(invA12),1,nfram).*R12);
    vsind = .5*ARA1 + repmat(invL12,1,nfram).*ARA2;
    vbind = repmat(invL22,1,nfram).*ARA2;
    
    % Enforcing the nonnegativity constraint (on vs only)
    neg = (vsind < 0) | (vbind < 0);
    vsind(neg) = 0;
    vs(:,:,ind) = vsind;
    vb(:,:,ind) = vbind;
end
spec = vs./vb;
     
return;



function spec = dsw_spec(hatRxx, f, d, tau_grid)

% DSW_SPEC Computes the SNR in all directions using the DS beamformer and
% frequency weighting
%
% spec = dsw_spec(hatRxx, f, d, tau_grid)
%
% Inputs:
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% d: microphone spacing in meters
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values
        
[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);
c = 343;
SINC = sinc(2*f*d/c);

SNR = zeros(nbin,nfram,ngrid);
for ind=1:ngrid,
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    SNR(:,:,ind) = repmat(-(1+SINC)/2,1,nfram) + repmat((1-SINC)/2,1,nfram).*(TR + 2*real(R12.*EXP))./(TR - 2*real(R12.*EXP));
end
spec = SNR;

return;



function spec = mvdrw_spec(hatRxx, f, d, tau_grid)

% MVDRW_SPEC Computes the SNR in all directions using the MVDR beamformer
% and frequency weighting
%
% spec = mvdrw_spec(hatRxx, f, d, tau_grid)
%
% Inputs:
% hatRxx : nbin x nfram x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% d: microphone spacing in meters
% tau_grid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nfram x ngrid array of SNR values
        
[nbin,nfram] = size(hatRxx(:,:,1,1));
ngrid = length(tau_grid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);
c = 343;
SINC = sinc(2*f*d/c);

SNR = zeros(nbin,nfram,ngrid);
for ind=1:length(tau_grid),
    EXP = repmat(exp(-2*1i*pi*tau_grid(ind)*f),1,nfram);
    NUM = real(R11.*R22 - R12.*R21)./(TR - 2*real(R12.*EXP));
    SNR(:,:,ind) = repmat(-(1+SINC)/2,1,nfram) + repmat((1-SINC)/2,1,nfram).*NUM./(.5*TR-NUM);
end
spec = SNR;

return;
