clear all
close all
Nt  = 50000;
ww = [];
for N =  [2500];
    for NrOfTrials = [1000];
        nK = 20;

        FPE_val = zeros(1, nK);
        AIC_val = zeros(1, nK);
        MDL_val = zeros(1, nK);

        FPE_order = zeros(NrOfTrials, 1);
        AIC_order = zeros(NrOfTrials, 1);
        MDL_order = zeros(NrOfTrials, 1);

        % The AR model
        omega = pi*[0.90 0.70 0.50 0.30 0.10];
        rho   =    [0.75 0.95 0.85 0.80 0.90];
        npairs = length(omega);
        roots_vec = zeros(2*npairs,1);
        roots_vec(1:npairs) = rho.*exp(1i*omega);
        roots_vec(npairs+1:end) = conj(roots_vec(1:npairs));
        AR = poly(roots_vec);

        for i=1:NrOfTrials

            % One realization of AR process of length Nt samples
            e = randn(Nt,1);
            y = filter(1, AR, e);

            % Pick last N samples from this AR process realization
            y_1 = y((end-N+1):end);
            y_1 = y_1+randn(N,1);
            FPE = [];
            AIC = [];
            MDL = [];
            % Order estimation loop
            for k_estim=1:nK

                % We use Yule-Walker method to obtain an estimate of the all-pole
                % AR parameters
                [~, E] = aryule(y_1,k_estim);

                FPE(k_estim) = ((N + k_estim)/(N-k_estim))*E
                AIC(k_estim) = N*log2(E) + 2*k_estim
                MDL(k_estim) = N*log2(E) + k_estim*log2(N)
                % Implement FPE, AIC, and MDL order selection criteria.
                % Use the variance E, obtained from aryule(), as the residual power
                % to drive each of the order selection criteria.
                % Store the result of each criteria in FPE_val(k_estim), AIC_val(k_estim), or
                % MDL_val(k_estim).

            end

            [~,FPE_order(i)] = min(FPE);
            [~,AIC_order(i)] = min(AIC);
            [~,MDL_order(i)] = min(MDL);
            % For each criteria, pick the value 'k*' which minimizes the criterion
            % and store the value of 'k*' in the corresponding vectors FPE_order(i),
            % AIC_order(i), or MDL_order(i).

        end

        % From vectors FPE_order, AIC_order, MDL_order, obtain the histogram of the
        % values of 'k*' and pick the one with the highest frequency as the order
        % estimate for the given criteria.

        hh_FPE = hist(FPE_order,1:nK);

        hh_AIC = hist(AIC_order,1:nK);

        hh_MDL = hist(MDL_order,1:nK);

        figure(1),clf,bar(1:nK,hh_FPE), title('FPE')
        figure(2),bar(1:nK,hh_AIC), title('AIC')
        figure(3),bar(1:nK,hh_MDL), title('MDL')

        [~,FPE_order_b] = max(hh_FPE)
        [~,AIC_order_b] = max(hh_AIC)
        [~,MDL_order_b] = max(hh_MDL)
        ww = [ww; N NrOfTrials  FPE_order_b AIC_order_b MDL_order_b]
    end
end
ww




