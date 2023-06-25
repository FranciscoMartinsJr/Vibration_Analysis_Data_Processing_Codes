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


%%
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



%% Harmonics detector
% clear input;
% clc;
% 
% userInput = input('Digite a frequência fundamental do sistema [Hz]: ', 's');
% userFrequency = str2double(regexp(userInput, '[\d.]+', 'match'));
% clear input;
% userInput2 = input('Digite a quantidade de harmonicas: ', 's');
% userHarmonics = str2double(regexp(userInput2, '[\d.]+', 'match'));
% 
% closestIndex = find(fmag >= userFrequency, 1, 'first');
% windowSize = 7500; % Size of the frequency window to consider
% windowIndices = max(closestIndex - windowSize, 1) : min(closestIndex + ...
%     windowSize, numel(fmag));
% [~, maxAmplitudeIndex] = max(MAG(windowIndices));
% 
% closestIndex = windowIndices(maxAmplitudeIndex);
% closestFrequency = fmag(closestIndex);
% closestAmplitude = MAG(closestIndex);
% 
% numHarmonics = userHarmonics;
% 
% detectedHarmonics = zeros(numHarmonics, 1);
% detectedAmplitudes = zeros(numHarmonics, 1);
% detectedFrequencies = zeros(numHarmonics, 1);
% 
% for i = 1:numHarmonics
%     harmonicIndex = closestIndex * i; % Index of the harmonic
%     detectedHarmonics(i) = i; % Store the harmonic number
%     detectedAmplitudes(i) = MAG(harmonicIndex); % Retrieve the amplitude
%     detectedFrequencies(i) = fmag(harmonicIndex); % Retrieve the frequency
% end
% 
% 
% 
% msg = sprintf('Harmonicas Detectadas Apontando Desbalanceamento Mecânico:\n');
% for i = 1:numel(detectedHarmonics)
%     msg = sprintf('%sHarmonica %d: Amplitude = %.5f m/s², Frequency = %.5f Hz\n', ...
%         msg, detectedHarmonics(i), detectedAmplitudes(i), detectedFrequencies(i));
% end
% msgbox(msg, 'Harmonicas Detectadas');


%% Harmonics Detector CORRECTED
clear input;
clc;

userInput = input('Digite a frequência fundamental do sistema [Hz]: ', 's');
userFrequency = str2double(regexp(userInput, '[\d.]+', 'match'));
% clear input;
% userInput2 = input('Digite a quantidade de harmonicas: ', 's');
% userHarmonics = str2double(regexp(userInput2, '[\d.]+', 'match'));

% Find the index of the closest frequency to the user input
closestIndex = find(fmag >= userFrequency, 1, 'first');

fftData = MAG;
frequencyAxis = fmag;

% Find the index of the highest amplitude in the vicinity of the closest frequency
windowSize = 7500; % Size of the frequency window to consider
windowIndices = max(closestIndex - windowSize, 1) : min(closestIndex + windowSize, numel(frequencyAxis));
[~, maxAmplitudeIndex] = max(fftData(windowIndices));

% Calculate the actual closest frequency and its corresponding amplitude
closestIndex = windowIndices(maxAmplitudeIndex);
closestFrequency = frequencyAxis(closestIndex);
closestAmplitude = fftData(closestIndex);

% Determine the number of harmonics you want to detect
numHarmonics = 33;
detectedHarmonics = zeros(numHarmonics, 1);
detectedAmplitudes = zeros(numHarmonics, 1);
detectedFrequencies = zeros(numHarmonics, 1);

% Calculate the harmonics based on the new closest frequency
for i = 1:numHarmonics
    harmonicFrequency = closestFrequency * i;
    
    % Find the index of the closest frequency to the harmonic
    harmonicIndex = find(frequencyAxis >= harmonicFrequency, 1, 'first');
    
    % Find the index of the highest amplitude in the vicinity of the harmonic frequency
    harmonicWindowIndices = max(harmonicIndex - windowSize, 1) : min(harmonicIndex + windowSize, numel(frequencyAxis));
    [~, harmonicMaxAmplitudeIndex] = max(fftData(harmonicWindowIndices));
    
    % Calculate the actual harmonic frequency and its corresponding amplitude
    harmonicIndex = harmonicWindowIndices(harmonicMaxAmplitudeIndex);
    detectedHarmonics(i) = i;
    detectedAmplitudes(i) = fftData(harmonicIndex);
    detectedFrequencies(i) = frequencyAxis(harmonicIndex);
end

% Filter out any detected frequencies that have lower amplitudes than the closest frequency
validIndices = detectedAmplitudes >= closestAmplitude;
detectedHarmonics = detectedHarmonics(validIndices);
detectedAmplitudes = detectedAmplitudes(validIndices);
detectedFrequencies = detectedFrequencies(validIndices);

% Display the detected harmonics, their amplitudes, and frequencies using a message box
msg = sprintf('Harmonicas Detectadas Apontando Irregularidade Mecânica:\n');
for i = 1:numel(detectedHarmonics)
    msg = sprintf('%sHarmonica %d: Amplitude = %.6f m/s², Frequency = %.4f Hz\n', ...
        msg, detectedHarmonics(i), detectedAmplitudes(i), detectedFrequencies(i));
end
msgbox(msg, 'Detected Harmonics');