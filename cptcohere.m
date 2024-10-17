function [CptC,varargout] = cptcohere(x,y,varargin)
%CPTCOHERE  Calculates all component coherence (CptC) values (ordinary
%component and cross component coherence) for a given multi-input,
%single-output system. 
% 
% Note: if any inputs are too highly correlated (e.g., two identical 
% inputs) or if any inputs have zero power at a frequency of interest then 
% the matrix to be inverted in estimation of the system's optimum frequency
% response functions will be singular and component coherence results will 
% be unreliable. 
%   CptC = CPTCOHERE(X,Y), returns the component coherence decomposition of
%   the system with input set X and output Y using Welch's averaged,
%   modified periodogram method. Component coherence decomposes output
%   power into components attributable to inputs directly (ordinary
%   component coherence) vs  to interference between input pairs (cross
%   component coherence). X is an m x q array where m is the number of
%   points in each input signal and q is the number of inputs to the
%   system. Y is an m x 1 vector. CptC is an n x q x q array, where n is 
%   the number of FFT points used to calculate the needed PSD estimates.
%   CptC(f,i,j) is the component coherence between inputs i & j and the 
%   output at frequency f (i=j is ordinary component coherence, i/=j is 
%   cross component coherence). 
%  
%   By default, X and Y are divided into eight sections with 50% overlap,
%   each section is windowed with a Hamming window, and eight modified
%   periodograms are computed and averaged. These settings can be modified
%   by calling the function with additional inputs as follows, where 
%   variable names and definitions are described in "help cpsd". 
%   
%   CptC = CPTCOHERE(X,Y,WINDOW)
%   CptC = CPTCOHERE(X,Y,WINDOW,NOVERLAP)
%   CptC = CPTCOHERE(X,Y,WINDOW,NOVERLAP,NFFT)
%   [CptC,w] = CPTCOHERE(___)
%   [CptC,f] = CPTCOHERE(___,Fs)
%   [CptC,w] = CPTCOHERE(X,Y,WINDOW,NOVERLAP,W)
%   [CptC,f] = CPTCOHERE(X,Y,WINDOW,NOVERLAP,F,Fs)
%
%   [CptC,__,H] = CPTCOHERE(X,Y,__), also returns an n x q array of the
%   optimum frequency response functions used to estimate component
%   coherence, where H(f,i) is the optimum frequency response function from
%   input i to the output at frequency f.
%
% Written by Dr. Steven Charles and Nolan Howes, 2024.
    
    % Manage additional inputs
    nargoutchk(1,3)
    narginchk(2,6)
    window = []; noverlap = []; f = []; fs = [];
    if nargin > 2
        window = varargin{1};
        if nargin > 3; noverlap = varargin{2}; end
        if nargin > 4; f = varargin{3}; end
        if nargin > 5; fs = varargin{4}; end
    end
    
    % Calculate needed PSDs
    Giy = cpsd(x,y,window,noverlap,f,fs); %input-output cpsd
    Gyy = cpsd(y,y,window,noverlap,f,fs); %output psd
    [Gij,f] = cpsd(x,x,window,noverlap,f,fs,'mimo'); %input spectral density matrix

    % Calculate optimum frequency response functions and component
    % coherence at each f
    H = zeros(size(Giy)); 
    CptC = zeros(size(Gij));
    for i = 1:length(f)
        Gijf = squeeze(Gij(i,:,:)); %Get Gij @ f
        [U,S,V] = svd(Gijf); %SVD for matrix inversion
        Hf = (V*(S\U')*Giy(i,:).'); %Optimum FRF @ f
        CptC(i,:,:) = (diag(conj(Hf))*Gijf*diag(Hf))/Gyy(i); %Component Coherence @ f
        H(i,:) = Hf.'; %Append opt. FRF to vector
    end

    % Manage additional outputs
    if nargout == 2
        varargout = {f};
    elseif nargout == 3
        varargout = {f,H};
    end
end
