clear all
close all


%% Choose an AR model with which to generate data

n= 8; % number of zeros of the AR polynomial

j=sqrt(-1);
r1= 0.5; om1=3*pi/4; % modulus and phase of the first two roots
r2 = 0.8; om2 = 3*pi/5; % modulus and phase of the next two roots
r3 = 0.9; om3 = pi/6;  % modulus and phase of the last two roots
r4 = 0.998; om4 = pi/12;  % modulus and phase of the last two roots
rootsi =[r1*exp(j*om1) r1*exp(-j*om1) r2*exp(j*om2) r2*exp(-j*om2) ...
    r3*exp(j*om3) r3*exp(-j*om3)  r4*exp(j*om4) r4*exp(-j*om4)]; % all roots
A = poly(rootsi); % the AR polynomial having the roots given in rootsi
nA = length(A);

N0 = 1024; % number of frequency points on the frequency axis (in range [0; 2*pi))
omega0 = 2*pi/N0; % the resolution in frequency
% the set of frequency points omega = ii*omega0
om = 0:omega0:(2*pi);

%% Compute the true spectrum of the AR process
%% Evaluate P(omega) = 1/|A(exp(i*omega))|^2

H = zeros(length(om),1);
for ii = 1:length(om)
    % H(ii) = exp(-j*om(ii)).^(0:(nA-1))*A';
    H(ii) = 0;
    omega = om(ii);
    for iAR = 1:nA
        H(ii) = H(ii) + A(iAR)*(exp(-j*omega))^(iAR-1);
    end
    P(ii) = 1 ./ abs(H(ii)).^2;
end

figure(1),clf
subplot(222)
plot(om(1:N0/2),log10(P(1:N0/2)))
grid
title('True power spectrum of AR process','Fontsize',14)
xlabel('frequency $\omega$','interpreter','Latex','Fontsize',12)
ylabel('$\log_{10} P(\omega)$','interpreter','Latex','Fontsize',12)

%% Questions:
% Q1: The highest power spectrum value P(omega) is at the angular frequency:
% a) omega = pi/6; b) omega = pi/ 3; c) omega = 3*pi/4; d) omega = pi/12.

%% Generate one realization of the AR process

N = 200000;  % Number of data samples to be generated
y = rand(n,1); % initialize the memory of the AR process randomly
e = randn(N,1);
for ii = (n+1):N
    y(ii) = -sum( y(ii-(1:n))' .* A(2:end)) +e(ii);
end
% Keep only the last (N0 = 1024) data samples (the rest are discarded to be
% sure that the transient part of the AR process ended)
y1 = y((1:N0) + end-N0); % this is the realization further processed for getting the spectrum estimate
% Compute the FFT of the N0 = 1024 data points
Y1 = fft(y1);
P1 = abs(Y1).^2/N0; % this is the periodogram estimate
om1 = om(1:(N0/2));

subplot(223)
plot(om1,log10(P1(1:N0/2)))
grid
title('Estimated spectrum from one realization','Fontsize',18)
xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)
ylabel('$\log \hat P(\omega)$','interpreter','Latex','Fontsize',16)

% Q2: The highest estimated power spectrum value P1(omega) is at
% an angular frequency close to
% a) omega = pi/6; b) omega = pi/12; c) omega = pi/3; d) omega = 3*pi/4.

subplot(224)
plot(om1,log10(P1(1:N0/2)),'b')
hold on
plot(om(1:N0/2),log10(P(1:N0/2)),'r-','LineWidth',3)
grid
title('Estimated overlapped with true spectrum','Fontsize',18)
xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)
ylabel('$\log \hat P(\omega)$, $\log \hat P(\omega)$','interpreter','Latex','Fontsize',16)


subplot(221)
NN = N0;
plot(y1(1:NN))
grid
title('One realization of AR process','Fontsize',18)
xlabel('  time $t$','interpreter','Latex','Fontsize',16)
ylabel('$y_t$','interpreter','Latex','Fontsize',16)

%% Questions:
% Q3: The plot of y1 in subplot(221) shows an underlying repetitive process,
% with a period T of about
% a) 12 samples; b) 30 samples;d) 24 samples; c) 3 samples

% Q4: The period you notticed at Q3 corresponds to the angular frequency having the
% a) lowest power; b) highest power c) DC term (i.e., omega = 0);


%% Generate multiple realizations of the AR proces, find the associated
%% FFT,and average over all realizations at each frequency value

N_realiz = 100;
NG0 = 1024;   
    om0 = 2*pi/NG0;
    Pmean = zeros(NG0,1);
    tic
    for i_realiz = 1:N_realiz
        % Simulate the AR process
        y = zeros(N,1);
        y(1:n)=rand(n,1);
        N= 200000;
        e= randn(N,1);
        for ii = (n+1):N
            y(ii) = -sum( y(ii-(1:n))' .* A(2:end)) +e(ii);
        end
        % Select NG0 points at the end of the AR data
        y1 = y((1:NG0) + end-NG0);
        % Compute the fast Fourier transform

        Y1 = fft(y1);
        P1 = abs(Y1).^2/NG0;

        Pmean = Pmean + P1;
    end
    toc
    Pmean = Pmean/N_realiz;
    om1 = om0*(1:NG0/2);


    figure(3),clf,
    plot(om1,log10(P1(1:NG0/2)),'-b','Linewidth',0.5),hold on
    plot(om1,log10(Pmean(1:NG0/2)),'-r','Linewidth',1),hold on
    grid
    title('Estimated spectrum over many realizations','Fontsize',14)
    xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)
    ylabel('$\log_{10} E[\hat{P}(\omega)], \log_{10} \hat{P_1}(\omega)$','interpreter','Latex','Fontsize',16)
    plot(om(1:N0/2),log10(P(1:N0/2)),'-k','Linewidth',2)

    figure(10),clf
    plot(om1,log10(Pmean(1:NG0/2)),'-r'),hold on
    plot(om(1:N0/2),log10(P(1:N0/2)),'-b')
    grid
    title(['Average of estimated spectrum Nrealiz =  ' num2str(N_realiz) '  Ndata = ' num2str(NG0)],'Fontsize',12)
    xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)
    ylabel('$\log_{10} E[\hat{P}(\omega)], \log_{10} P(\omega)$','interpreter','Latex','Fontsize',16)


%Q5: By generating independent realizations and averaging the
% spectrum estimate obtained at a given omega, as we do in Pmean,
% one gets a method similar to
% a) Bartlett method b) Welch method c) Daniell method


%% Generate multiple realizations of the AR proces, find the associated
%% FFT,and average over all realizations at each frequency value

N_realizi = [100 100 400 400]; % Numbers of realizations to be used in experiments
NG0i =[1024 4096 1024 4096]; % Numbers of data to be used in the experiment
for ij = 1:4
    tic
    N_realiz = N_realizi(ij);
    NG0 = NG0i(ij);
    om0 = 2*pi/NG0;
    Pmean = zeros(NG0,1);
    tic
    for i_realiz = 1:N_realiz
        % Simulate the AR process
        y = zeros(N,1);
        y(1:n)=rand(n,1);
        N= 200000;
        e= randn(N,1);
        for ii = (n+1):N
            y(ii) = -sum( y(ii-(1:n))' .* A(2:end)) +e(ii);
        end
        % Select NG0 points at the end of the AR data
        y1 = y((1:NG0) + end-NG0);
        % Compute the fast Fourier transform

        Y1 = fft(y1);
        P1 = abs(Y1).^2/NG0;
        Pmean = Pmean + P1;
    end
    toc
    Pmean = Pmean/N_realiz;
    om1 = om0*(1:NG0/2);


    figure(10+ij),clf
    plot(om1,log10(Pmean(1:NG0/2)),'-r'),hold on
    plot(om(1:N0/2),log10(P(1:N0/2)),'-b')
    grid
    title(['Average of estimated spectrum  Nrealiz =  ' num2str(N_realiz) '  Ndata = ' num2str(NG0)],'Fontsize',12)
    xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)
    ylabel('$\log_{10} E[\hat{P}(\omega)], \log_{10} P(\omega)$','interpreter','Latex','Fontsize',16)

    % Find corresponding points on the omega axis
    Pmean_alligned = zeros(N0/2,1);
    P_alligned = zeros(N0/2,1);
    for ik = 1:(N0/2)
        [~, kstar] = min( abs(om1-om(ik)) );
        Pmean_alligned(ik) = Pmean(kstar);
        P_alligned(ik) = P(ik);
    end

    figure(20+ij),clf
    %figure,clf
    plot(om(1:N0/2),log10(Pmean_alligned)-log10(P_alligned),'-r'),hold on
    grid
    title(['Difference Nrealiz =  ' num2str(N_realiz) '  Ndata = ' num2str(NG0)],'Fontsize',12)
    xlabel('frequency $\omega$','interpreter','Latex','Fontsize',16)

end



%Q6: The spectrum estimate Pmean is very close to the ideal spectrum P,
% however, there are differences between the two. How one can ensure that
% P and Pmean become closer:
% a) Use a larger number of realizations N_realiz; b) Use a larger value of NG0;
% c) Use a smaller value of NG0; d) Use a smaller number of realizations
% N_realiz;e) Increasing N_realiz produces closer P to Pmean, than when
% increasing NG0

