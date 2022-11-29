clear all
close all
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

%% Task 1. Run RLS for the signals from Room 1 and Room 2

%% Calling RLS using Matlab's DSP toolbox


M = 200;         % filter length
lambda = 1; % forgetting factor, use values in the range [0.998,1]

% Recursive least squares from Matlab's DSP toolbox
hRLS = dsp.RLSFilter(...
    'Length',M,...
    'Method','Conventional RLS',...
    'ForgettingFactor',lambda);
%% Task 1.1 % Run the RLS for the signals in Room 1
% Run the RLS for the signals in Room 1
[y,e] = step( hRLS,v,s1 );

% Compute the MSE between the cleaned speech e(n) and the clean speech c(n)
Ns =length(s1);
interval= (round(Ns/4):Ns);
MSE_RLS_room1 = sum((c(interval)-e(interval)).^2)/length(interval)


%% Task 1.2 % Run the RLS for the signals in Room 2
% Complete in the lines bellow the missing code 
[y2,e2] = step(hRLS,v,s2);
Ns2 =length(s2);
interval2= (round(Ns2/4):Ns2);
MSE_RLS_room2 = sum((c(interval2)-e2(interval2)).^2)/length(interval2)

%% Task 2. Try RLS with 8 equally space values of lambda between 0.9995 and 1.
% Complete in the lines bellow the missing code 

MSE_RLS_room1_s = [];
lambda_s = linspace(0.9995,1,8);
for i1 = 1:length(lambda_s)
    lambda = lambda_s(i1);
    % Recursive least squares from Matlab's DSP toolbox
    hRLS = dsp.RLSFilter(...
        'Length',M,...
        'Method','Conventional RLS',...
        'ForgettingFactor',lambda);
    [y,e] = step( hRLS,v,s1 );
    MSE_RLS_room1_s(i1) = sum((c(interval)-e(interval)).^2)/length(interval)
end

MSE_RLS_room2_s = [];
lambda_s = linspace(0.9995,1,8);
for i1 = 1:length(lambda_s)
    lambda = lambda_s(i1);
    % Recursive least squares from Matlab's DSP toolbox
    hRLS = dsp.RLSFilter(...
        'Length',M,...
        'Method','Conventional RLS',...
        'ForgettingFactor',lambda);
    [y,e] = step( hRLS,v,s2 );
    MSE_RLS_room2_s(i1) = sum((c(interval)-e(interval)).^2)/length(interval)
end

figure(1),clf, plot(lambda_s, MSE_RLS_room1_s,'or-')
hold on
plot(lambda_s, MSE_RLS_room2_s,'vb-')
grid on
xlabel('lambda'),ylabel('MSE RLS')
legend('Room 1','Room 2')

%% Q4. For which lambda do we obtain the lowest MSE when s(t) equals speech and noise through room 1.wav
%% Q5. For which lambda do we obtain the lowest MSE when s(t) equals speech and noise through room 2.wav

format long
lambda_s

%% Task 3 Least Squares: Run the file LS.m and load matrix A


% You are given the Matlab ﬁle LS.m which loads the matrix A and the necessary audio signals.
% Run the file LS.m
LS
figure(3),imagesc(A(1:400,:)),colormap(gray),colorbar
figure(4),imagesc(A(31600+(1:400),:)),colormap(gray),colorbar

%% Q6. Inspect the structure of A. Is the matrix built using pre-windowing, post-windowing method, or pre- and post-windowing?

%% Task 4 Least Squares: Run the file LS.m to get initial data
%% Implement the LS estimator given in Equation 7  (in the pdf Exercise5.pdf)

% Consider the case of Room 1
d = s1;
%% Q7. % Compare the following estimators
w1 = (A'*A)\(A'*d);
w2 = A\d;
w3 = inv(A'*A)*(A'*d);
% First compare the maximum absolute difference
max( abs(w1-w2))
max( abs(w1-w3))
% And then compare the relative absolute difference
max( abs(w1-w2))/max(abs(w1))
max( abs(w1-w3))/max(abs(w1))

% Q7. Which of the following statements are true?
% - w1 = (A'*A)\(A'*d) gives very similar results to w2 = A\d;
% - w3 = inv(A'*A)*(A'*d) is the ideal solution but is very different from w1 = (A'*A)\(A'*d) 
% - The solution w2 = A\d; involves less numerical errors when solving the
% system Aw=d in the least square sense

%% Task 5. Compute the cleaned speech using the LS parameters
%% Implement the LS filtering given in Equation 7 and 8 (in the pdf Exercise5.pdf)

% The case of Room 1
d = s1;
w2 = A\d;
e = d-A*w2;
interval= 1:Ns;
MSE_LS_room1= sum((c(interval)-e(interval)).^2)/length(interval)

% The case of Room 2
d = s2;
w2 = A\d;
e = d-A*w2;
interval= 1:Ns;
MSE_LS_room2= sum((c(interval)-e(interval)).^2)/length(interval)

%% Task 6. Least squares noise cancellation with segmentation

% Follow the instruction in the pdf and write your own solution

% The case of Room 1, LS for the first segment
M = 200;
d = s1;
A1 = [];
A2 = [];
segment = 2;
N = length(v)/segment;
Ns = length(s1)/segment;
d1 = d(1:1:N)
d2 = d(N+1:end)
    %for ind = length(d)/segment
   for ii = 1:M
       A1 = [A1 [zeros(ii-1,1); v(1:(N-ii+1))]];
   end

    %w1 = A_temp\d(1:1:NN,1);
    w1 = (A1'*A1)\(A1'*d1);
    %w1 = A1\d1;
    e1 = d1-A1*w1;
    %w2 = (A1'*A1)\(A1'*d1);
   
   for ii = 1:M
       A2 = [A2 [zeros(ii-1,1); v(N+1:(2*N-ii+1))]];
   end

    %w2 = A2\d2;
    w2 = (A2'*A2)\(A2'*d2);
    e2 = d2-A2*w2;
    e = [e1;e2];
    %end
Ns =length(s1);
%MSE_RLS_room1 = sum((c(interval)-e(interval)).^2)/length(interval)
%MSE_RLS_room1 = sum((e-c).^2)/length(c)
% after that, answer to Questions Q9,10,11,12
interval= 1:Ns;
MSE_LS_room2= sum((c(interval)-e(interval)).^2)/length(interval)






