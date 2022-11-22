%% Exercise 1 Random signals
% We study the distribution of the values of audio signals. First the basic MATLAB functions are introduced for generating
% data according to the uniform and normal distribution. Secondly we look at audio data and try to describe it in probabilistic
% terms. Your task is to run the provided code and understand its functionality and modify it when asked. Return your answer
% using the course Moodle page.

%% 1 a: MATLAB generated data
% In MATLAB the most used random data generators are the rand and randn functions. A short description of the functions is
% given below:
%
% * rand(n, m) generates a _n×m_ matrix containing pseudo random values drawn from the standard uniform dis-
% tribution on the open interval (0, 1);
% * randn(n, m) returns a _n×m_ matrix containing pseudo random values drawn from the standard normal distri-
% bution with zero mean and unit variance.
%
% Run the two m-functions listed below and answer to the questions marked
% in the Moodle Quize Exercise 1 a:
%
% # uniform_distribution.m
% # normal_distribution.m

%% Exercise 1 b: Real audio data
% For this experiment we use a fragment of the audio data BeeGees.wav. The goal is to see if the data has a known
% distribution; to see if a probabilistic description can be found. In this
% case we try to ﬁt a normal distribution on the data and see how well the
% ideal distribution agrees with the data.
% The description and code for the exercise is in the matlab function real_random_signals.m  

