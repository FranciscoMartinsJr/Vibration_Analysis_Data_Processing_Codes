%% Vibration Analysis Data Processing Script

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

x = x_wdc - mean(x_wdc);
y = y_wdc - mean(y_wdc);
z = z_wdc - mean(z_wdc);

%% Convert x, y, and z data from 'g' to 'm/s²'

g = 9.80665;
x = x.*g;
y = y.*g;
z = z.*g;


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

plot_samples = ((0:length(x_wdc)-1)./(60*fs))';

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


%% Combining all axis

mag = (x + y + z);

Nmag = length(mag);
MAG = (abs(fft(mag, Nmag)))./Nmag; % Normalized FFT
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


%% Amplitudes Detector 
clear input;
clc;

% userInput = input('Digite a frequência fundamental do sistema [Hz]: ', 's');
% userFrequency = str2double(regexp(userInput, '[\d.]+', 'match'));
% clear input;
% userInput2 = input('Digite a quantidade de harmonicas: ', 's');
% userHarmonics = str2double(regexp(userInput2, '[\d.]+', 'match'));
userFrequency = 15;

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
numHarmonics = 3;
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
msg = sprintf('Amplitudes Detectadas Apontando Irregularidade Mecânica:\n');
for i = 1:numel(detectedHarmonics)
    msg = sprintf('%sAmplitude = %.6f m/s², Frequencia = %.4f Hz\n', ...
        msg, detectedAmplitudes(i), detectedFrequencies(i));
end
msgbox(msg, 'Amplitudes Detectadas');   

