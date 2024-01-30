function SINR_output = F_SINR(L, K, H, H_hat, rho_ul, pilots)

SINR_output = zeros(K, L);
for l = 1:L
    A = H_hat(:,:,l)';

    for k = 1:K
        ICI1 = 0;        % Intra-cell interference                                               
        for k1 = 1:K
            if k1~=k
                if pilots(:,k,l) == pilots(:,k1,l)
                    ICI1 = ICI1 + rho_ul*norm(A(k,:)*H(:,k1,l,l))^2;
                end
            end
        end

        ICI2 = 0;         % Inter-cell interference                                              
        for j1 = 1:L
            for k1 = 1:K
                if j1~=l 
                    if pilots(:,k,l)==pilots(:,k1,j1)
                        ICI2 = ICI2 + rho_ul*norm(A(k,:)*H(:,k1,l,j1))^2;
                    end
                end
            end
        end

        ICIN = norm(A(k,:))^2;     % equivalent noise                                     
        SINR_output(k, l) = rho_ul*norm(A(k,:)*H(:,k,l,l))^2/(ICI1 + ICI2 + ICIN);

    end
end

end