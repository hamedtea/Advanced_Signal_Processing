% clean speech
[c,fs] = audioread('clear_speech.wav', 'native');
c = double(c);

% unpertubed noise
[v,fs] = audioread('noise_source1.wav', 'native');
v = double(v);

% room 1
[s1,fs] = audioread('speech_and_noise_through_room_11.wav', 'native');
s1 = double(s1);

% room 2
[s2, fs] = audioread('speech_and_noise_through_room_22.wav', 'native');
s2 = double(s2);

% filter length
M = 200;
N = length(v);

% Construct the data matrix A
A = [];
for ii = 1:M
    A = [A [zeros(ii-1,1); v(1:(N-ii+1))]];
end
save('A.mat',"A")
% load data matrix A
AA = load('A.mat')
A = AA.A;