%% Uniform distribution.

% Generate a vector u1 containing 100000 values from a pseudorandom variable generated using a uniform
% distribution.
u1 = rand(100000, 1);

% The function hist can be used to split the range of the values into M=50
% bins and then count the number of values that fall within each of the
% defined bins
M = 50;

figure(1),clf 
histogram(u1, M);
ylabel('Bin count'); 
xlabel('Values');
sum( (u1>0.06) & (u1<0.08))
%% Uniform distribution in matab
pd1 = makedist('Uniform');
x = 0:0.01:1;
y = pdf(pd1, x);
hist(x,y)
%% Q.1a1: Which of the following statements are true regarding the plot in Figure 1? (multiple choices might be correct)
%
% * The 4th bin in Figure 1 corresponds to the interval with edges 0.06 and
% 0.08 and the number of values of u1 falling in the bin is close to
% 100000/50, agreeing with the assumption of uniform probability distribution
% * The number of values that fall into the bin 4 is not equal to
% 100000/50, hence violating the assumption of uniform distribution of the
% matlab function rand

%------------------------------------------------------------------------
% Generate B, as a 100x1000 matrix containing 100000 pseudorandom variables from an
% uniform distribution.
B = rand(100, 1000);

% We use the same hist function to count the values of the first row from B.
figure(2),clf
histogram(B(:, :), M);
ylabel('Bin count'); 
xlabel('Values');

%% Q.1a2: Which of the following statements are true regarding the plot in Figure 2? (multiple choices might be correct)
%
% * All values in the first row of B about are contained in the interval (0,1)
% * The dispersion of the number of values in each of the 50 bins is very
% high and we cannot make any statement about the distribution used for
% generating them, based only on the plotted histogram 
% * By using all elements of the matrix B for drawing a histogram with 50
% bins, one can be more reassured about the uniform distribution of the
% values into the bins.

%------------------------------------------------------------------------
% Generate the vector u2  containing 100000 pseudorandom values from an uniform
% distribution U(a,b)
u2 = -4 + 6 * rand(100000000, 1);

% The function hist is used to count the number of values.
figure(3),clf
hist(u2, M);
ylabel('Bin count'); 
xlabel('Values');
sum(u2 > 3)

%% Q.1a3: Which of the following statements are true regarding the plot in Figure 3? (multiple choices might be correct)
%
% * The uniform distribution has parameters a=-4 and b=6
% * The uniform distribution has parameters a=-4 and b = 2
% * By generating more than 100000 values, there is a  high chance to obtain at 
%  least one value larger than 3
 

%% 
pd = makedist('Normal','mu',75,'sigma',10);
x = 25:125;
y = pdf(pd, x);
figure
plot(x,y)
grid
