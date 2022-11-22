% Advanced Signal Processing
% Exercise 2 Optimal Wiener filters
%

clear
close all

% load the audio sample that is treated as the input signal
[d, fs] = audioread('Tamara_Laurel_-_Sweet_extract.wav', 'native');
d = double(d);

% simulate the channel output with infinite signal-to-noise ratio
[u, w_true] = simulate_channel(d, Inf);

% extract a segment from both signals
s_start = 8;
d = d(s_start*fs+1:(s_start+1)*fs);
u = u(s_start*fs+1:(s_start+1)*fs);

%% Use xcorr() to obtain the biased and unbiased autocorrelation estimates
% of u. Use maximum lag value 300.
tau_max = 300;
%[p,lags] = autocorr(u,NumLags=tau_max-1); %econometrics add-in
[p,lags] = xcorr(u, tau_max, 'unbiased'); 
%p_unbiased = p/(length(u)-tau_max);
figure
%plot(0:1:tau_max-1,p(299:1:end),'ro')
stem(lags,p)
figure
plot(lags,p)

%% Write your Wiener filter implementation here using d and u. 
% Use unbiased estimate for cross-correlation vector p (xcorr(d,u)) 
% and matrix R, both with maximum lag value 5. Use toeplitz() for
% reordering the autocorrelation estimate into the matrix R.
% Compare the obtained w with the channel parameters w_true.
ta_max = 5;
[p,lags] = xcorr(d,u, ta_max, 'unbiased');
[rr,~] = xcorr(u, ta_max, 'unbiased'); 
r = toeplitz(rr(ta_max+1:end));
figure
stem(lags,p);
figure
stem(lags,rr);
%ww=inv(r)*p(ta_max+1:end);
w = r\p(ta_max+1:end);
d_hat = filter(w,1,u);
e = sum((d-d_hat).^2)/length(u); %mean square error













