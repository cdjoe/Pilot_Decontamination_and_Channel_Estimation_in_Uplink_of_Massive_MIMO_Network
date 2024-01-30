function h_hat = CE_aid(M, K, L, H, pilots, rho_pilot)

h_hat = zeros(M,K,L);   % least square estimator
% for i = 1:L
%     for k = 1:K
%         h_hat(:,k,i) = (rx_pilot(:,:,i)*pilots(:,k,i)) / (sqrt(rho)*pilot_N);
% %        h_hat(:,k,i) = rx_pilot(:,:,i)*conj(pilots(:,k)) / (pilots(:,k).' * conj(pilots(:,k)));
%     end
% end


    for i = 1:L
        for k = 1:K
            h_hat(:,k,i) = H(:,k,i,i);                                  % effective channel signal
            
            for i1 = 1:L
                for k1 = 1:K
                    
                    if pilots(:,k1,i1) == pilots(:,k,i) 
                        h_hat(:,k,i) = h_hat(:,k,i) + H(:,k1,i,i1);     % inter-cell interference of channel estimation
                    end
                    
                end
            end
            
            h_hat(:,k,i) = h_hat(:,k,i) + 1/sqrt(rho_pilot)*random('norm', 0, 1, M, 1);
        end
    end
        
end