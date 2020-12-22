function [X] = tremor(dat, y, Fs, Fmax, n, vb )
%% Function
% Tremor detection
% Detect tremor based on range and RMS amplitude signal
%
%% USAGE
%   tremor(dat, y, Fs, Fmax, n, vb);
%
%% INPUTS
%   dat (1xn struct)	MAVIS data structure containing EMA data
%   y (string)          name of EMA trace to analyze
%   Fs (int)            sampling frequency
%   Fmax (tin)          the maximum frequency of the DFT
%   n (int)             the number of points used to calculate the fft
%   vb (int)            optionally plot FFT (vb == 1)


%% EXAMPLE USAGE : tremor(trial01, "TB", 400, 20, 2048, 0); 

%% Created:   22-jun-16	Michael Proctor (mike.i.proctor@gmail.com)
%  Modified:  01-sep-20    Jidde Jacobi

    names = {dat.NAME};
    ix1 = strcmp( names, y);
    sensheight = dat(ix1).SIGNAL(:,3);  % select z-axis of EMA trace of interest
    
    hpf = designfilt('highpassiir', 'FilterOrder',8, 'PassbandFrequency',1, 'PassbandRipple',0.05, 'SampleRate', 400); %create HPF
    filt = filtfilt( hpf, sensheight);  % apply HPF


    if (nargin < 4)
        n	= 512;                      % no points FFT
    end

    m = n/2;
    nqst        = Fs/2;					% calculate Nyquist frequency
    X1          = fft(filt,n);			% calculate DFT on signal y
    Pxx         = X1 .* conj(X1)/n;		% power density spectrum
    Pxx1        = 10*log10(Pxx);        % Db scale
    fbase       = Fs*(1:m)/n;			% create frequency axis
    nmax        = m*Fmax/nqst;
    rootms      = rms(filt);
    [pks,locs]  = findpeaks(Pxx(1:round(nmax)), fbase(1:round(nmax)), 'MinPeakProminence',0.75*(max(Pxx)-min(Pxx)));   % retrieve locations of peaks

     
    if vb
        figure; plot(fbase(1:round(nmax)),Pxx(1:round(nmax)),locs,pks,'or'); 
        title("Spectrum")
        xlabel('Frequency (Hz)');
        ylabel("Intensity");
        xlim([0 8])
    end

    
    X.filename          = dat; 
    X.pklocs            = locs;
    X.rms               = rootms;
    X.fourier           = X1;
    X.power_density     = Pxx;
    X.power_density_dB  = Pxx1;
    X.frequency         = fbase;

    
   
end