function h_hat = CE_LS(M, K, L, rx_pilot, pilots, rho, pilot_N)

h_hat = zeros(M,K,L);   % least square estimator
for i = 1:L
    for k = 1:K
        h_hat(:,k,i) = (rx_pilot(:,:,i)*pilots(:,k,i)) / (sqrt(rho)*pilot_N);
%        h_hat(:,k,i) = rx_pilot(:,:,i)*conj(pilots(:,k)) / (pilots(:,k).' * conj(pilots(:,k)));
    end
end
        
end