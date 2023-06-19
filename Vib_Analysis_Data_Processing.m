%% Conversion from txt to csv

data = readmatrix('LavRapida_NvH20Alto_PesoNormal.txt');


%% Matrix N:1 samples

samples = ((0:length(data)-1)./1000)';

%% Append sample value to cvs data

data = [samples data];


%% Slice cvs data to obtain t, x, y, and variables


t = data(:,1); % Valores do tempo
x = data(:,2); % Valores das acelerações eixo x
y = data(:,3); % Valores das acelerações eixo y
z = data(:,4); % Valores das acelerações eixo z

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


%% FFTs calculation (in m/s² and dB)

clc;

Over_Sample = 1;

N1 = Over_Sample.*length(x);

X1 = abs(fft(x, N1))./N1;
magX1 = 20*log10(abs(fft(x, N1)));

Y1 = abs(fft(y, N1))./N1;
magY1 = 20*log10(abs(fft(y, N1)));

Z1 = abs(fft(z, N1))./N1;
magZ1 = 20*log10(abs(fft(z, N1)));

f1 = (0: fs/N1 : fs - fs/N1)';

%% FFTs plots in m/s² 

close all; clc;

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

%% Transpose

x = x';
y = y';
z = z';

%% Samples interpolation to slice data in specific places

t1 = 2257540;
t2 = 2484570;

t_s = (t1:t2)';
x_s = interp1(x, t_s);
y_s = interp1(y, t_s);
z_s = interp1(z, t_s);


%% Sliced data plots in respect to time in minutes

close all;

sliced_samples = ((0:length(x_s)-1)./(60*fs))';

figure
subplot(2,2,1)
plot(sliced_samples, x_s, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo x');

subplot(2,2,2)
plot(sliced_samples, y_s, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo y');

subplot(2,2,3)
plot(sliced_samples, z_s, '-k')
xlabel('Tempo [min]');
ylabel('Amplitude [m/s²]');
title('Dados do Acelerômetro para o eixo z');

%% FFTs calculation for sliced data

N_s = length(x_s);

Over_Sample_s = 1;

N_s = Over_Sample_s.*N_s;

X1_s = abs(fft(x_s, N_s))./N_s;
magX1_s = 20*log10(abs(X1_s));

Y1_s = abs(fft(y_s, N_s))./N_s;
magY1_s = 20*log10(abs(Y1_s));

Z1_s = abs(fft(z_s, N_s))./N_s;
magZ1_s = 20*log10(abs(Z1_s));

f1_s = (0 : fs/N_s : fs - fs/N_s)';

%% FFT plots in m/s²
clc; close all;
y_max = 100;

figure
subplot(2,2,1)
plot(f1_s, X1_s, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo x');

subplot(2,2,2)
plot(f1_s, Y1_s, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo y');

subplot(2,2,3)
plot(f1_s, Z1_s, '-k')
xlim([0 fn]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [m/s²]');
title('Espectro de Fourier para o eixo z');

%% FFT plots in dB
clc; close all;
figure
subplot(3,1,1)
plot(f1_s, magX1_s, '-k')
xlim([0 fn]);
ylim([0 y_max]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo x');

subplot(3,1,2)
plot(f1_s, magY1_s, '-k')
xlim([0 fn]);
ylim([0 y_max]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo y');

subplot(3,1,3)
plot(f1_s, magZ1_s, '-k')
xlim([0 fn]);
ylim([0 y_max]);
xlabel('Frequência [Hz]');
ylabel('Amplitude [dB]');
title('Espectro de Fourier para o eixo z');
