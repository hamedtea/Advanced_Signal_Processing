clear
close all

% load the audio sample that is treated as the input signal
[d, fs] = audioread('Tamara_Laurel_-_Sweet_extract.wav', 'native');
d = double(d);

% simulate the channel output
[u, w_true] = simulate_channel(d, 80);

% extract a segment from both signals
s_start = 8;
d = d(s_start*fs+1:(s_start+1)*fs);
u = u(s_start*fs+1:(s_start+1)*fs);

% estimate autocorrelation and cross-correlation
N_THETA = 6;
r = xcorr(u, u, N_THETA-1, 'unbiased');
r = r(N_THETA:end);
R = toeplitz(r);

p = xcorr(d, u, N_THETA-1, 'unbiased');
p = p(N_THETA:end);

% solve the Wiener filter
w_wiener = inv(R)*p;

w = zeros(N_THETA,1);

% use steepest descent to find the wiener filter

%mus = [1e-10 1e-9 1e-8 1e-7 3e-7]; %  different step sizes. 
% the last one is close to the largest stable step size.
mus = [1e-6];
figure(1),clf

for mu = mus
    
    N = 5000; % number of iterations
    
    % Use steepest descent. 
    % Columns of Wt are the filter coefficients at different time instants (or iterations).
    Wt = steepest_descent(R, p, w, mu, N);
        
    % compute the mean square error in the filter coefficients
    % at different time instants compared to the wiener filter
    mse_coeffs = zeros(N,1);
    for n = 1:N
        mse_coeffs(n) = mean((Wt(:,n)-w_wiener).^2); 
    end
    plot(mse_coeffs), grid on, hold on
    
end
xlabel('Iteration')
ylabel('MSE in the coefficients')
legend(num2str(mus'))

%% Q5.  What is the value of the mean square error after 100 iterations with µ = 10−7 ?
%% Q6. Does increasing the step size parameter µ increase or decrease the convergence speed?
% (Try other µ values as well)

%% Q7. What happens if you use a step size of µ = 10−6 ?















