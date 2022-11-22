function w_save = steepest_descent(R, p, w, mu, N)
% The steepest descent algorithm.
%
% IN:
%   R - the autocorrelation matrix
%   p - the cross-correlation vector
%   w - the initial filter coefficients
%   mu - the step size
%   N - maximum number of iterations
%
% OUT:
%   w_save - the filters at different iterations as columns(!) of w_save
%            i.e. the final filter is at w_save(:,end)

w_save = zeros( numel(p), N+1 );
w_save(:,1) = w;
for t = 1:N
    
     w = w  + mu*(p-R*w);% [implement the filter update here];
     w_save(:,t+1) = w;

end