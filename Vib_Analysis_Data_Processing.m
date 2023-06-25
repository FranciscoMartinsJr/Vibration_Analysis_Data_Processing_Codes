%% Vibrational Analysis Data Processing Script

%% Conversion from txt to csv

data = readmatrix('LavRapida_NvH20Alto_PesoNormal.txt');

%% Matrix N:1 samples

samples = ((0:length(data)-1)./1000)';

%% Append sample value to cvs data

data = [samples data];


%% Slice cvs data to obtain t, x, y, and variables


t = data(:,1); % Valores do tempo
x_wdc = data(:,2); % Valores das acelerações eixo x
y_wdc = data(:,3); % Valores das acelerações eixo y
z_wdc = data(:,4); % Valores das acelerações eixo z

%% Sample time and frequency

Ts = mean(diff(t)); %Periodo de amostragem médio em segundos
fs = 1/Ts;          %Frequência de amostragem
fn = fs/2;          %Nyquist

%% Remove DC values from x, y and z data

x = x - mean(x);
y = y - mean(y);
z = z - mean(z);

%% Convert x, y, and z data from 'g' to 'm/s²'

g = 9.80665;
x = x.*g;
y = y.*g;
z = z.*g;


%% Slice undesired/random collected data from the beggining

% NOTE: The first high acceleration amplitude is from the first
%       abrupt break when the washing machine turns to settle
%       the clothes inside so it shouldn't be claimed as noise.

t1 = 12238;
t2 = 2502096;

t_s = (t1:t2)';
x = interp1(x, t_s);
y = interp1(y, t_s);
z = interp1(z, t_s);

%% PLot x, y, and z data values in respect to sample

close all;

plot_samples = (0:length(x)-1)';

figure
subplot(3,1,1)
plot(plot_samples, x, '-k')
xlabel('Amostras');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo x');

subplot(3,1,2)
plot(plot_samples, y, '-k')
xlabel('Amostras');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo y');

subplot(3,1,3)
plot(plot_samples, z, '-k')
xlabel('Amostras');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo z');

%% Plot x, y and z data values in respect to time in minutes

close all;

plot_samples = ((0:length(x)-1)./(60*fs))';

figure
subplot(3,1,1)
plot(plot_samples, x, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo x');

subplot(3,1,2)
plot(plot_samples, y, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo y');

subplot(3,1,3)
plot(plot_samples, z, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo z');

%% Plot in time in the same figure
close all;

plot_samples = ((0:length(x)-1)./(60*fs))';

figure
plot(plot_samples, x_wdc)
hold on
plot(plot_samples, y_wdc)
hold on
plot(plot_samples, z_wdc)
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo z');
legend('Eixo x', 'Eixo y', 'Eixo z');

%% FFTs calculation (in m/s² and dB)

clc;

N1 = length(x);

X1 = (abs(fft(x, N1)))./N1; % Normalized FFT
imagX1 = (imag(fft(x, N1)))./N1;
magX1 = 20*log10(abs(fft(x, N1)));

Y1 = (abs(fft(y, N1)))./N1;
imagY1 = (imag(fft(y, N1)))./N1;
magY1 = 20*log10(abs(fft(y, N1)));

Z1 = (abs(fft(z, N1)))./N1;
imagZ1 = (imag(fft(z, N1)))./N1;
magZ1 = 20*log10(abs(fft(z, N1)));

f1 = (0: fs/N1 : fs - fs/N1)';

%% FFTs plots in m/s²
clc;

y_lim = 80;

figure
subplot(2,2,1)
plot(f1, X1, '-k')
xlim([0 fn]);
yticks('auto');
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo x');

subplot(2,2,2)
plot(f1, Y1, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo y');

subplot(2,2,3)
plot(f1, Z1, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo z');

figure
subplot(2,2,1)
plot(f1, imagX1, '-k')
xlim([0 fn]);
yticks('auto');
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo x');

subplot(2,2,2)
plot(f1, imagY1, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo y');

subplot(2,2,3)
plot(f1, imagZ1, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo z');

%% FFT plot in dB
close all; clc;

figure
plot(f1, magX1, '-k')
xlim([0 fn]);
ylim([0 y_lim]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo x');

figure
plot(f1, magY1, '-k')
xlim([0 fn]);
ylim([0 y_lim]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo y');

figure
plot(f1, magZ1, '-k')
xlim([0 fn]);
ylim([0 y_lim]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo z');

%% Combining all axis

mag = (x + y + z);
Nmag = length(mag);
MAG = (abs(fft(mag, Nmag)))./Nmag; % Normalized FFT
MAG_IMAG = (imag(fft(mag, Nmag)))./Nmag;
powerMAG = 20*log10((abs(fft(mag, Nmag)))./Nmag);
fmag = (0: fs/Nmag : fs - fs/Nmag)';

mag_vs_time = ((0:Nmag-1)./(60*fs))';
mag_samples = (0:Nmag-1)';


%% Ploting in frequency

clc; close all

figure
plot(fmag, MAG, '-k')
xlim([0 fn]);
% set(gca, 'YScale', 'log')
yticks('auto');
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier');

figure
plot(fmag, MAG_IMAG, '-k')
xlim([0 fn]);
% set(gca, 'YScale', 'log')
yticks('auto');
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier');

%% Plotting in time

clc; close all;

figure
plot(mag_vs_time, mag, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo x');

%% Warning string when certain amplitude exceeded the threshold
clc

f_range1 = 0;
f_range2 = 500;
amplitude = 0.005; % Amplitude threshold
% Index the frequency range
frequencyIndices = find(fmag >= f_range1 ...
    & fmag <= f_range2);

% Get amplitude values from indexed frequencies
amplitudeValues = MAG(frequencyIndices, :);

% Get the exceeding frequencies in regards to amplitude threshold
exceedingIndices = find(amplitudeValues > amplitude);

% Select from indexed freqs the frequencies
selectedFrequencies = f1(frequencyIndices);

% Get which frequencies met the amplitude threshold
exceedingFrequencies = selectedFrequencies(exceedingIndices);

% Get the list associated with the exceeding frequencies
exceedingAmplitudes = amplitudeValues(exceedingIndices);

% Select the maximum amplitude between the list
[~, maxAmplitudeIndex] = max(exceedingAmplitudes);
highestAmplitude = exceedingAmplitudes(maxAmplitudeIndex);

% Select at which frequency is associated with that amplitude
correspondingFrequency = exceedingFrequencies(maxAmplitudeIndex);

% Create a string representation of the variables
varStr = sprintf(['Desbalanceamento mecânico detectado.\n' ...
    'Amplitude [m/s²]: %f\nFrequência [Hz]: %f'], highestAmplitude, ...
    correspondingFrequency);

% Display the variables using msgbox
msgbox(varStr, 'Atenção');

% %% Low-pass Filter
% clc;
% % Filter parameters
% cutoff_freq = 1;  % Cutoff frequency for the low-pass filter (Hz)
% order = 4;  % Order of the Butterworth filter
% 
% % Design the Butterworth filter
% nyquist_freq = fs/2;  % Nyquist frequency
% normalized_cutoff = cutoff_freq/nyquist_freq;  % Normalize cutoff frequency
% [b, a] = butter(order, normalized_cutoff, 'low');  % Design the Butterworth filter coefficients
% 
% % Apply the filter to the signal
% x_filtered = filter(b, a, x);
% y_filtered = filter(b, a, y);
% z_filtered = filter(b, a, z);
% 
% %% FFTs calculation of filtered signal (in m/s² and dB)
% 
% clc;
% 
% Nfiltered = length(x_filtered);
% 
% XFILTERED = (abs(fft(x_filtered, Nfiltered)))./Nfiltered; % Normalized FFT
% imagXFILTERED = (imag(fft(x_filtered, Nfiltered)))./Nfiltered;
% magXFILTERED = 20*log10(abs(fft(x_filtered, Nfiltered)));
% 
% YFILTERED = (abs(fft(y_filtered, Nfiltered)))./Nfiltered;
% imagYFILTERED = (imag(fft(y_filtered, Nfiltered)))./Nfiltered;
% magYFILTERED = 20*log10(abs(fft(y_filtered, Nfiltered)));
% 
% ZFILTERED = (abs(fft(y_filtered, Nfiltered)))./Nfiltered;
% imagZFILTERED = (imag(fft(y_filtered, Nfiltered)))./Nfiltered;
% magZFILTERED = 20*log10(abs(fft(y_filtered, Nfiltered)));
% 
% fFILTERED = (0: fs/Nfiltered : fs - fs/Nfiltered)';

% %% FFTs low-pass filtered plots in m/s²
% clc; close all;
% 
% y_lim = 80;
% 
% figure
% subplot(2,2,1)
% plot(f1, XFILTERED, '-k')
% xlim([0 fn]);
% yticks('auto');
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo x');
% 
% subplot(2,2,2)
% plot(f1, YFILTERED, '-k')
% xlim([0 fn]);
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo y');
% 
% subplot(2,2,3)
% plot(f1, ZFILTERED, '-k')
% xlim([0 fn]);
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo z');
% 
% figure
% subplot(2,2,1)
% plot(f1, imagXFILTERED, '-k')
% xlim([0 fn]);
% yticks('auto');
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo x');
% 
% subplot(2,2,2)
% plot(f1, imagYFILTERED, '-k')
% xlim([0 fn]);
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo y');
% 
% subplot(2,2,3)
% plot(f1, imagZFILTERED, '-k')
% xlim([0 fn]);
% xlabel('Frequência [Hz]');
% ylabel('Amplitude [m/s²]');
% title('Espectro de Fourier para o eixo z');