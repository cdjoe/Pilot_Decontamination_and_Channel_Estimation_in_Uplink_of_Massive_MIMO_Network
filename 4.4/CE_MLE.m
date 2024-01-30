function h_hat = CE_MLE(M, K, L, rx_pilot, pilots, rho, pilot_N, Beta)

LS = zeros(M,K,L);   % least square estimator
for i = 1:L
    for k = 1:K
        LS(:,k,i) = (rx_pilot(:,:,i)*pilots(:,k,i)) / (sqrt(rho)*pilot_N);
%        LS(:,k,i) = Y(:,:,i)*conj(S(:,k)) / (S(:,k).' * conj(S(:,k)));
    end
end

zeta_mle = zeros(L,K);
for i = 1:L
    for k = 1:K
        zeta_mle(i,k) = (norm(LS(:,k,i))^2) / M;
    end
end
        
h_hat = zeros(M,K,L);
for i = 1:L
    for k = 1:K
        h_hat(:,k,i) = M *Beta(k,i,i) /(norm(LS(:,k,i))^2) *LS(:,k,i);
    end
end
        
end