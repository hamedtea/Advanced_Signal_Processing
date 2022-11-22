%% Behind the scenes of Scenario 1
% function s1_sim = BehindScenesScenario1

% % Check the values in Wt_s1(k,i) which was used to simulate the signal s1
% 
% load the filter Wt_s1 accounting for wall reflections in room 1
load('Wt.mat');

% clear (clean) speech
[c,fs] = audioread('clear_speech.wav', 'native');
c = double(c);
N = length(c);

% corrupting noise
v = 700*randn(N,1);
% room 1 response to the source of noise, induced at microphone
% Wt_s1(:,1)'
corrupting_noise1 = filter(Wt_s1(:,1),1,v);

% recorded signal at microphone
s1_sim = c+corrupting_noise1;
player = audioplayer(s1_sim/max(abs(s1_sim)),fs);
playblocking(player);

[max(corrupting_noise1) min(corrupting_noise1) max(s1_sim) min(s1_sim) max(c) min(c)]

% visualize the variation of Wt_s1 along time dimension 
% The order of the true filter (used in simulation)
K = size(Wt_s1,1);
N = size(Wt_s1,2);
% The true weights at the initial time (Wt_s1(1:K,1))
% The true weights at the final time Wt_s1(:,N)
figure(1),clf, plot(1:K, Wt_s1(1:K,1),'or-', 1:K, Wt_s1(1:K,N),'vb-')
xlabel('k, index of the weight')
title({'The initial weights of the simulation filter and','final weights of the simulation filter'})

% Are the weights changing along the trajectory in timw?
Max_True_Weights = max(Wt_s1,[],2);
Min_True_Weights = min(Wt_s1,[],2);

figure(2),clf, plot(1:K, Min_True_Weights,'or-', 1:K, Max_True_Weights,'vb-')
xlabel('k, index of the weight')
title({'The minimum of each weight along time of the simulation and',...
    'The maximum of each weight along time of the simulation'});
 
if(0)
    audiowrite('noise_source1.wav',int16(v),fs)
    audiowrite('speech_and_noise_through_room_11.wav',int16(s1_sim),fs)
end
