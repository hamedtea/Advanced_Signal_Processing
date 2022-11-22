%% Behind the scenes of Scenario 2
% function BehindScenesScenario2
[c,fs] = audioread('clear_speech.wav', 'native');
c = double(c);
[max(c) min(c)]
N= length(c);
% initial source of noise
[v,fs] = audioread('noise_source1.wav', 'native');
v = double(v);
% % Check the values in Wt_s2(k,i) which was used to simulate the signal s1
% 
load('Wt.mat');
Wt_s2(:,1)'
% at time n use the filter coefficients in Wt_s1(:,n)
vect_v= zeros(M,1); % initialize the input vector with zeros
corrupting_noise2 = zeros(N,1);
for n = 1:N
    vect_v = [v(n); vect_v(1:(M-1))]; % create the vector [u(n); u(n-1); ... ;u(n-M)]
    corrupting_noise2(n) = vect_v'*Wt_s2(:,n);
end
[max(corrupting_noise2) min(corrupting_noise2) max(s2) min(s2)]


s2_sim = c+corrupting_noise2;
player = audioplayer(s2_sim/max(abs(s2_sim)),fs);
playblocking(player);

% The order of the true filter (used in simulation)
K = size(Wt_s2,1);
N = size(Wt_s2,2);
% The true weights at the initial time (Wt_s1(1:K,1))
% The true weights at the final time Wt_s1(:,N)
figure(1),clf, plot(1:K, Wt_s2(1:K,1),'or-', 1:K, Wt_s2(1:K,N),'vb-')
xlabel('k, index of the weight')
legend('initial','final')
title({'The initial weights of the simulation filter and','final weights of the simulation filter'})

% Are the weights changing along the trajectory in time?
Max_True_Weights = max(Wt_s2,[],2);
Min_True_Weights = min(Wt_s2,[],2);

figure(2),clf, plot(1:K, Min_True_Weights,'or-', 1:K, Max_True_Weights,'vb-')
xlabel('k, index of the weight')
legend('minimum','maximum')
title({'The minimum of each weight along time of the simulation and',...
    'The maximum of each weight along time of the simulation'});

if(0)
    audiowrite('speech_and_noise_through_room_22.wav',int16(s2_sim),fs)
end

