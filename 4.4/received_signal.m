function Y = received_signal(M, N, L, rho, H, pilots)

Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
% rho = 10^(q/10);   % the uplink power / TX SNR(dB)
sumG = zeros(M,N,L);

for i = 1:L
    for l = 1:L
        sumG(:,:,i) = sumG(:,:,i) + H(:,:,i,l)*pilots(:,:,i)';
    end
    Y(:,:,i) = Y(:,:,i) + sqrt(rho)*sumG(:,:,i) + (randn(M,N)+1i*randn(M,N))/sqrt(2);
%             Y(:,:,i) = Y(:,:,i) + sqrt(q)*sumG(:,:,i) ;
end
        
end