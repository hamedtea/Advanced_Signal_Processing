%% Real audio signal
clear;
clc;
close all;
play_sounds = 1;

[Y1,FS] = audioread('BeeGees.wav', 'native');

Y = double(Y1);

% Play in audioplayer the waveform in Y.
if(play_sounds)
    player = audioplayer(Y1(1:15*FS),FS);
    playblocking(player);
end

% The number of elements in signal Y, i.e., signal length.
nY = numel(Y);

% We clean zero valued samples from the signal to make the fitting more
% illustrative. You can also try with this line commented to see what happens.
Y(Y==0) = [];


%% Fit a normal distribution to the signal

% Now use normfit() to compute the mean ('muhat') and variance ('sigmahat') of the values of Y:

[muhat,sigmahat] = normfit(Y);

fitting_precision = 2^5;

Xboundaries = (-2^15):fitting_precision:(2^15-1);

Xcenters =  Xboundaries + fitting_precision/2;

% Cumulative normal distribution at each boundary point P(x<Xb(j)):
P = normcdf(Xboundaries, muhat, sigmahat);

% Probability of each interval:
Pint = P(2:end) - P(1:(end-1));

% obtain histogram of the signal values
[hh,edge] = histcounts(Y,Xboundaries);

figure(1),clf

% Plot the histogram of the signal values
plot(Xcenters(2:end),hh),hold on
% Plot the histogram obtained from the fitting of the normal distribution
plot(Xcenters(2:end), nY*Pint, '--r');
ylabel('Value count');
xlabel('Centers values');

title('full audio signal')

legend('empirical histogram','fitting of normal distribution')

%% Select a small segment containing only one second of audio.

Tbase = 1000000;
Ya = Y(Tbase+(1:44100));

nYa = numel(Ya);

% Now use normfit() to compute the mean ('muhat') and variance ('sigmahat') of the values of Ya:

[muhat_a,sigmahat_a] = normfit(Ya);

fitting_precision = 2^5;

Xboundaries = (-2^15):fitting_precision:(2^15-1);

Xcenters =  Xboundaries + fitting_precision/2;

% Cumulative normal distribution at each boundary point P(x<Xb(j)):
P = normcdf(Xboundaries, muhat_a, sigmahat_a);

% Probability of each interval:
Pint = P(2:end) - P(1:(end-1));

% obtain histogram of the signal values
[hh,edges] = histcounts(Ya,Xboundaries);

figure(2),clf

% Plot the histogram of the signal values
plot(Xcenters(2:end),hh),hold on
% Plot the histogram obtained from the fitting of the normal distribution
plot(Xcenters(2:end), nYa*Pint, '--r');
ylabel('Value count');
xlabel('Centers values');

title('one second of audio')

legend('empirical histogram','fitting of normal distribution')

%% Signal quatization

% Uniform quantizer with step size Delta = 2^9.
Delta = 2^9;
% Assume the centers to be:
Xcenters = -2^15:2^9:2^15;

% The quantizer is uniform, the reconstructed value is
% round(Y(i)/Delta)*Delta
% For example 1:1 the ramp signal -1000:10000 is quantized to :
true_y = -10000:10000;
quantized_y = round(true_y/Delta) * Delta;

figure(3),clf
plot(true_y,true_y,'-b'),hold on
plot(true_y,quantized_y,'--r')
grid on
legend('true signal','quantized signal')

% Let's quantize now the audio signal:
quantized_Y = round(Y/Delta)*Delta;
quant_error = Y-quantized_Y;

figure(4),clf
hist(quant_error, min(quant_error):max(quant_error));


