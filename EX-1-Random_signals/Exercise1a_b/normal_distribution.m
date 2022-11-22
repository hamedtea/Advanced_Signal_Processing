%% Normal distribution.

% Generate a vector n1 containing 100000 pseudorandom values from a normal
% distribution; the call is similar with that of rand:
n1 = randn(100000, 1);


% The function hist can be used to count the number of values that fit in
% M=50 bins.
M = 60;

figure(1),clf
hist(n1, M);
ylabel('Value count'); 
xlabel('Values');
[mean(n1) var(n1)]
%% Q.1a4: Which of the following statements are true regarding the plot in Figure 1?
%
% * The mean of the distribution is 0
% * The standard deviation of the distribution is 1
% * The variance of the distribution is 3
 


%---------------------------------------
% Generate the vector n2 using randn function in the following way:
m = 4.23;
s = 3.123;
n2 = m + s * randn(100000, 1);

figure(2),clf
hist(n2, M);
ylabel('Value count');
xlabel('Values');
[mean(n2) sqrt(var(n2))]

%% Q.1a5: Which of the following statements are true regarding the plot in Figure 2?
%
% * The mean of the distribution is 4.23
% * The standard deviation of the distribution is 3.123
% * The variance of the distribution is 20

