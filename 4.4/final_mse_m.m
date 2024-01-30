% Channel estimation for massive MIMO TDD systems assuming pilot contamination and flat fading
% Rayleigh fading

clear all;
% close all;
%% Parameter
L = 7;   % cells
K = 10;   % users

ant_M = 10:10:200;   % antenna elements
Mm = length(ant_M);

R_Cell = 1000;
r_Min = 100;
alpha = 3.8;
sigma_shadow = 8;

a = 0.05;   % the cross-cell interference level
lambda = 0.1;
% q = 10;

n = 300;   % iterations
mean_eta_mle = zeros(4,Mm,n);
mean_eta_ls = zeros(4,Mm,n);
mean_eta_mmse = zeros(4,Mm,n);

final_eta_mle = zeros(4,Mm);
final_eta_ls = zeros(4,Mm);
final_eta_mmse = zeros(4,Mm);

for num = 1:n
    display(num);
    
    for ii = 1:Mm
        M = ant_M(ii);
%         display(M);
        
        %% Channel coefficients      
%         H = zeros(M,K,L,L);   % channel coefficients
%         Beta = zeros(K,L,L);   % large-scale coefficient
% 
%         for i = 1:L
%             for l = 1:L
%                 for k = 1:K
%                     if i == l
%                         Beta(k,i,i) = 1;
%                     else
%                         Beta(k,i,l) = a;
%                     end
%                     H(:,k,i,l) = H(:,k,i,l) + sqrt(Beta(k,i,l))*(randn(M,1)+1i*randn(M,1))/sqrt(2);   % *small-scale coefficient CN(0,1)
%                 end
%             end
%         end
        
        [H, Beta] = F_H_Generate(M, L, K, R_Cell, r_Min, sigma_shadow, alpha);
        
        %% Generate Pilots (Zadoff-Chu sequence)
        S = 10;
        N = S;   % the length of the pilot sequence (N>=K)
        temp = zeros(L,S);
        for i=1:L
            temp(i,:) = randperm(S);
        end
        
        pilot = F_ZC(N,S);        
        pilots = zeros(N,K,L);
        for l=1:L
            for k=1:K
                pilots(:,k,l) = pilot(:,temp(i,k));
            end
        end
            
        %--- Received Pilot Signal ---
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        q = 10^(10/10);   % the uplink power / TX SNR(dB)
        Y = received_signal(M,N,L,q,H,pilots);

        %--- MLE channel estimator ---
        MLE = CE_MLE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MLE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mle(1,ii,num) = sum(sum_m) / L;
        
        %--- LS channel estimator ---
        LS = CE_LS(M,K,L,Y,pilots,q,N);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(LS(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_ls(1,ii,num) = sum(sum_m) / L;
        
        %--- MMSE channel estimator ---
        MMSE = CE_ideal_MMSE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MMSE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mmse(1,ii,num) = sum(sum_m) / L;
        
        %% soft
        [pilots, N] = F_SPRS2(L,K,Beta,lambda);        
            
        % Received Pilot Signal
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        q = 10^(10/10);   % the uplink power / TX SNR(dB)
        Y = received_signal(M,N,L,q,H,pilots);
        
        %--- MLE channel estimator ---
        MLE = CE_MLE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MLE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mle(2,ii,num) = sum(sum_m) / L;
        
        %--- LS channel estimator ---
        LS = CE_LS(M,K,L,Y,pilots,q,N);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(LS(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_ls(2,ii,num) = sum(sum_m) / L;
        
        %--- MMSE channel estimator ---
        MMSE = CE_ideal_MMSE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MMSE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mmse(2,ii,num) = sum(sum_m) / L;
        
        %% wgcpa
        [P, N] = F_WGCPA(L,K,S,Beta);        
            
        % Received Pilot Signal
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        q = 10^(10/10);   % the uplink power / TX SNR(dB)
        Y = received_signal(M,N,L,q,H,P);

        %--- MLE channel estimator ---
        MLE = CE_MLE(M,K,L,Y,P,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MLE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mle(3,ii,num) = sum(sum_m) / L;
        
        %--- LS channel estimator ---
        LS = CE_LS(M,K,L,Y,P,q,N);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(LS(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_ls(3,ii,num) = sum(sum_m) / L;
        
        %--- MMSE channel estimator ---
        MMSE = CE_ideal_MMSE(M,K,L,Y,P,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MMSE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mmse(3,ii,num) = sum(sum_m) / L;
        
        %% soft+wgcpa
        [pilots, N] = F_SPRS_WGCPA2(L,K,Beta,lambda,S);
        
            
        % Received Pilot Signal
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        q = 10^(10/10);   % the uplink power / TX SNR(dB)
        Y = received_signal(M,N,L,q,H,pilots);

        %--- MLE channel estimator ---
        MLE = CE_MLE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MLE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mle(4,ii,num) = sum(sum_m) / L;
        
        %--- LS channel estimator ---
        LS = CE_LS(M,K,L,Y,pilots,q,N);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(LS(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_ls(4,ii,num) = sum(sum_m) / L;
        
        %--- MMSE channel estimator ---
        MMSE = CE_ideal_MMSE(M,K,L,Y,pilots,q,N,Beta);
        
        eta_MLE = zeros(L,K);
        sum_m = zeros(L,1); 
        for i = 1:L
            for k = 1:K  
                eta_MLE(i,k) = (norm(MMSE(:,k,i)-H(:,k,i,i))^2) / M;   % MSE per antenna of the MMSE estimator
            end
            sum_m(i,1) = sum(eta_MLE(i,:)) / K;
        end

        mean_eta_mmse(4,ii,num) = sum(sum_m) / L;
                
    end
        
end

for ii=1:Mm
    for num=1:n
        final_eta_mle(1,ii) = final_eta_mle(1,ii) + mean_eta_mle(1,ii,num)/n;
        final_eta_mle(2,ii) = final_eta_mle(2,ii) + mean_eta_mle(2,ii,num)/n;
        final_eta_mle(3,ii) = final_eta_mle(3,ii) + mean_eta_mle(3,ii,num)/n;
        final_eta_mle(4,ii) = final_eta_mle(4,ii) + mean_eta_mle(4,ii,num)/n;
        
        final_eta_ls(1,ii) = final_eta_ls(1,ii) + mean_eta_ls(1,ii,num)/n;
        final_eta_ls(2,ii) = final_eta_ls(2,ii) + mean_eta_ls(2,ii,num)/n;
        final_eta_ls(3,ii) = final_eta_ls(3,ii) + mean_eta_ls(3,ii,num)/n;
        final_eta_ls(4,ii) = final_eta_ls(4,ii) + mean_eta_ls(4,ii,num)/n;
        
        final_eta_mmse(1,ii) = final_eta_mmse(1,ii) + mean_eta_mmse(1,ii,num)/n;
        final_eta_mmse(2,ii) = final_eta_mmse(2,ii) + mean_eta_mmse(2,ii,num)/n;
        final_eta_mmse(3,ii) = final_eta_mmse(3,ii) + mean_eta_mmse(3,ii,num)/n;
        final_eta_mmse(4,ii) = final_eta_mmse(4,ii) + mean_eta_mmse(4,ii,num)/n;
    end
end

%% Plot
figure;

plot(ant_M,final_eta_ls(1,:),'c-x'); hold on
plot(ant_M,final_eta_ls(2,:),'b-x'); hold on
plot(ant_M,final_eta_ls(3,:),'g-x'); hold on
plot(ant_M,final_eta_ls(4,:),'r-x'); hold on

plot(ant_M,final_eta_mle(1,:),'c-o'); hold on
plot(ant_M,final_eta_mle(2,:),'b-o'); hold on
plot(ant_M,final_eta_mle(3,:),'g-o'); hold on
plot(ant_M,final_eta_mle(4,:),'r-o'); hold on

plot(ant_M,final_eta_mmse(1,:),'c-*'); hold on
plot(ant_M,final_eta_mmse(2,:),'b-*'); hold on
plot(ant_M,final_eta_mmse(3,:),'g-*'); hold on
plot(ant_M,final_eta_mmse(4,:),'r-*'); 

legend('Random:LS','SPRS:LS','WGC-PA:LS','SPRS+WGC-PA:LS','Random:MLE','SPRS:MLE','WGC-PA:MLE','SPRS+WGC-PA:MLE','Random:MMSE','SPRS:MMSE','WGC-PA:MMSE','SPRS+WGC-PA:MMSE');
axis([10 200 0.2 0.45]);
title('MSE versus M under varying pilot assignments');
xlabel('number of BS collocated antennas,M'); ylabel('MSE');
grid on;




































