function h_hat = CE_ideal_MMSE(M, K, L, rx_pilot, pilots, rho, pilot_N, Beta)

sum_beta = 0;
zeta = zeros(L,K);
for i = 1:L
    for k = 1:K
        for l = 1:L
            sum_beta = sum_beta + Beta(k,i,l);             
        end
        zeta(i,k) = sum_beta + (1/rho/pilot_N);
        sum_beta = 0;
    end
end

LS = zeros(M,K,L);   % least square estimator
for i = 1:L
    for k = 1:K
        LS(:,k,i) = (rx_pilot(:,:,i)*pilots(:,k,i)) / (sqrt(rho)*pilot_N);
%        LS(:,k,i) = Y(:,:,i)*conj(S(:,k)) / (S(:,k).' * conj(S(:,k)));
    end
end
        
h_hat = zeros(M,K,L);
for i = 1:L
    for k = 1:K
        h_hat(:,k,i) = LS(:,k,i) *Beta(k,i,i) /zeta(i,k);
    end
end
        
end