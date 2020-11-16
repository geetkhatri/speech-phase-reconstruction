% prompt the user to enter filename
filename = input('\n\nFilename: ', 's');
while exist(filename, 'file') ~= 2
    filename = input('File does not exist. Try again: ', 's');
end

info = audioinfo(filename);
if info.NumChannels ~= 1
    error('Audio must be mono.')
end

tic

% Input signal
[x, fs] = audioread(filename);

fsr = 8000;
%fsr = 8192;                         % frequency of resampled signal
M = 256;                            % window size
N = 256;                            % FFT size
H = M / 8;                          % hop size
w = hamming(M);                     % Hamming window of size M
a = 0.54;                           % constant = 0.54 for Hamming window

w = H * w / sum(w);                 % normalise window

% Resampling
if fs > fsr
    [D, L] = rat(fs / fsr);         % decimation and interpolation factors
    xr = resample(x, L, D);         % resampled signal
    fs = fsr;
end

% Add white Gaussian noise to signal at 20 dB SNR
xr = awgn(xr, 20);

% STFT
[X, f, t] = spectrogram(xr, w, M - H, N, fs);       % STFT of xr
magX = abs(X);                                      % magnitude spectrum
phaX = princ(angle(X));                             % phase spectrum

nframes = length(t);                % number of frames
nbands = length(f);                 % number of frequency bands
binwidth = f(2) - f(1);             % FFT bin width

% Fundamental frequencies
cd voicebox
% [f0, tv, pv, fv] = fxpefac(xr,fs,[0.0469 H/fs 10.5781],'g'); % fundamental frequencies
% f0 = [zeros(8, 1); f0; zeros(8, 1)];
% pv = [zeros(8, 1); pv; zeros(8, 1)];
p.fres = 1.81 * fs / M;
[f0, tv, pv, ~] = fxpefac(x, fs, H/fs, 'G', p);  % fundamental frequencies
voiced = zeros(nframes, 1);         % tells whether a frame is voiced
for i = 1 : nframes
    if pv(i) >= 0.9
        voiced(i) = 1;
    end
end
cd ../

% Phase reconstruction
phaY = phaX;               % reconstructed phase, initialized to phaX
for l = 2 : nframes
    
    % voiced sound
    if voiced(l) == 1
        
       % number of harmonics
       nharm = floor(fs ./ (2 * f0(l)));
       
       % harmonic frequencies
       fh = (1 : nharm).' * f0(l);
       
       % harmonic frequency closest to bin frequency, for each k
       fk = zeros(nbands, 1);      
       for k = 1 : nbands
           [magdif, I] = min(abs(fh - f(k)));
           fk(k) = fh(I);
       end
       
       % not onset of voiced sound
       if voiced(l - 1) ~= 0
           % normalised fk in radians
           wk = 2 * pi * fk / fs;
       
           % phase reconstruction along time
           phaY(:, l) = princ(phaY(:, l - 1) + wk * H);
       end
       
       % fk mapped to bin index notation
       bk = fh * N / fs;
       
       % bins that contain harmonic frequencies
       bh = round(bk);
       
       % extent of leakage of harmonic components into neighboring bands
       delk = ceil(fk(1) * N / (2 * fs));
       
       % phase reconstruction along frequency
       for i = [-delk : -1, 1 : delk]
           for j = 1 : nharm
               if 1 + bh(j) + i <= nbands
                   phaY(1 + bh(j) + i, l) = princ(phaY(1 + bh(j), l) ...
                          - phiW(2 * pi * (bh(j) - bk(j)) / N, a, M) ...
                          + phiW(2 * pi * (bh(j) - bk(j) + i) / N, a, M));
               end
           end
       end
       
    end
   
end

% Reconstructed STFT
Y = magX .* exp(1i * phaY);

% Reconstructed signal
[y, ty] = iosr.dsp.istft(Y, N, H, fs);      % inverse STFT
y = real(y);                                % neglect imaginary parts

toc

% Plot magnitude spectrum
figure
spectrogram(xr, w, M - H, N, fs, 'yaxis')
title('Magnitude of noisy and reconstructed signals')

% Plot phase spectrum of noisy signal
figure
surf(t, f / 1000, phaX, 'EdgeColor', 'None')
view(2)
title('Phase of noisy signal')
xlabel('Time (secs)')
ylabel('Frequency (kHz)')

% Plot phase of reconstructed signal
figure
surf(t, f / 1000, phaY, 'EdgeColor', 'None')
view(2)
title('Phase of reconstructed signal')
xlabel('Time (secs)')
ylabel('Frequency (kHz)')