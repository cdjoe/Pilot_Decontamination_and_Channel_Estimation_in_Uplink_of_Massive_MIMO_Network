clear;

clc; 
tic;
% system parameters
M = 512;    % number of BS antennas M
L = 7;     % cells
K = 10;     % users
S = 15;     % number of pilots

lambda=0.1;

R_Cell = 1000;
r_Min = 100;
alpha = 3.8;
sigma_shadow = 8;
% q = 15;
rho_ul = 10^(15/10);     % 15 dB                                                        

% Test_num = 300;
Test_num = 15;
UL_SINR_LS_sw = zeros(K, L, Test_num);
UL_SINR_MMSE_sw = zeros(K, L, Test_num);
UL_SINR_MLE_sw = zeros(K, L, Test_num);

UL_SINR_LS_wg = zeros(K, L, Test_num);
UL_SINR_MMSE_wg = zeros(K, L, Test_num);
UL_SINR_MLE_wg = zeros(K, L, Test_num);

UL_SINR_LS_sp = zeros(K, L, Test_num);
UL_SINR_MMSE_sp = zeros(K, L, Test_num);
UL_SINR_MLE_sp = zeros(K, L, Test_num);

UL_SINR_LS_ra = zeros(K, L, Test_num);
UL_SINR_MMSE_ra = zeros(K, L, Test_num);
UL_SINR_MLE_ra = zeros(K, L, Test_num);
SS = zeros(1,Test_num);

% Simulation
for i_test = 1:Test_num
        display(i_test);
                
        % Generate Channel vector
        [H, Beta] = F_H_Generate(M, L, K, R_Cell, r_Min, sigma_shadow, alpha);
        
        %% Generate Pilots (Zadoff-Chu sequence)
        N = S;   % the length of the pilot sequence (N>=K)
        pilot = F_ZC(N,S);
        
        %% random  
        temp = zeros(L,S);
        for i=1:L
            temp(i,:) = randperm(S);
        end
                
        pilots = zeros(N,K,L);
        for l=1:L
            for k=1:K
                pilots(:,k,l) = pilot(:,temp(i,k));
            end
        end
            
        %--- Received Pilot Signal ---
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        Y = received_signal(M,N,L,rho_ul,H,pilots);
        
        H_est = CE_aid(M, K, L, H, pilots, rho_ul);  
        
        UL_SINR_LS_ra(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        % ===== LS --> SINR =====%
%         H_est = CE_LS(M, K, L, Y, pilots, rho_ul, N);  
%         
%         UL_SINR_LS_ra(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
%                
%         % ===== MMSE --> SINR =====%
%         H_est = CE_ideal_MMSE(M, K, L, Y, pilots, rho_ul, N, Beta);  
%                 
%         UL_SINR_MMSE_ra(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        % ===== MLE --> SINR =====%
        H_est = CE_MLE(M, K, L, Y, pilots, rho_ul, N, Beta);  
               
        UL_SINR_MLE_ra(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        
        %% sprs
        [pilots, N] = F_SPRS2(L,K,Beta,lambda);
        
        Y = received_signal(M,N,L,rho_ul,H,pilots);
        
        H_est = CE_aid(M, K, L, H, pilots, rho_ul);  
        
        UL_SINR_LS_sp(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        % ===== LS --> SINR =====%
%         H_est = CE_LS(M, K, L, Y, pilots, rho_ul, N);  
%         
%         UL_SINR_LS_sp(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
%                
%         % ===== MMSE --> SINR =====%
%         H_est = CE_ideal_MMSE(M, K, L, Y, pilots, rho_ul, N, Beta);  
%                 
%         UL_SINR_MMSE_sp(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        % ===== MLE --> SINR =====%
        H_est = CE_MLE(M, K, L, Y, pilots, rho_ul, N, Beta);  
               
        UL_SINR_MLE_sp(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        %% wgcpa
        [P, N] = F_WGCPA(L,K,S,Beta);
        
        Y = received_signal(M,N,L,rho_ul,H,P);
        
        H_est = CE_aid(M, K, L, H, P, rho_ul);  
        
        UL_SINR_LS_wg(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,P);
        
        % ===== LS --> SINR =====%
%         H_est = CE_LS(M, K, L, Y, P, rho_ul, N);  
%         
%         UL_SINR_LS_wg(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,P);
%                
%         % ===== MMSE --> SINR =====%
%         H_est = CE_ideal_MMSE(M, K, L, Y, P, rho_ul, N, Beta);  
%                 
%         UL_SINR_MMSE_wg(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,P);
        
        % ===== MLE --> SINR =====%
        H_est = CE_MLE(M, K, L, Y, P, rho_ul, N, Beta);  
               
        UL_SINR_MLE_wg(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,P);
        
        %% SPRS+WGC-PA       
        [pilots, N] = F_SPRS_WGCPA2(L,K,Beta,lambda,S);
        
         %===== Received Pilot Signal =====%
%         Y = zeros(M,N,L);   % received uplink training sequences at i-th BS
        Y = received_signal(M,N,L,rho_ul,H,pilots);
        
        H_est = CE_aid(M, K, L, H, pilots, rho_ul);  
        
        UL_SINR_LS_sw(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
        
        % ===== LS --> SINR =====%
%         H_est = CE_LS(M, K, L, Y, pilots, rho_ul, N);  
%         
%         UL_SINR_LS_sw(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
%                
%         % ===== MMSE --> SINR =====%
%         H_est = CE_ideal_MMSE(M, K, L, Y, pilots, rho_ul, N, Beta);  
%                 
%         UL_SINR_MMSE_sw(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
%         
        % ===== MLE --> SINR =====%
        H_est = CE_MLE(M, K, L, Y, pilots, rho_ul, N, Beta);  
               
        UL_SINR_MLE_sw(:, :, i_test) = F_SINR(L,K,H,H_est,rho_ul,pilots);
            
        SS(i_test) = SS(i_test) + N;
end

Bw = 1;
mu0=0.05;
T_1=1-S/K*mu0;
% T_2=1-(SS+3)/K*mu0;



for i_test = 1:Test_num
    T_soft=1-SS(i_test)/K*mu0;
    
    for l=1:L
        for k=1:K
            UL_SINR_LS_ra(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_LS_ra(k,l,i_test));
            UL_SINR_LS_sp(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_LS_sp(k,l,i_test));
            UL_SINR_LS_wg(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_LS_wg(k,l,i_test));
            UL_SINR_LS_sw(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_LS_sw(k,l,i_test));
            
%             UL_SINR_MMSE_ra(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_MMSE_ra(k,l,i_test));
%             UL_SINR_MMSE_sp(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_MMSE_sp(k,l,i_test));
%             UL_SINR_MMSE_wg(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_MMSE_wg(k,l,i_test));
%             UL_SINR_MMSE_sw(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_MMSE_sw(k,l,i_test));
%             
            UL_SINR_MLE_ra(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_MLE_ra(k,l,i_test));
            UL_SINR_MLE_sp(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_MLE_sp(k,l,i_test));
            UL_SINR_MLE_wg(k,l,i_test)=Bw*T_1*log2(1+UL_SINR_MLE_wg(k,l,i_test));
            UL_SINR_MLE_sw(k,l,i_test)=Bw*T_soft*log2(1+UL_SINR_MLE_sw(k,l,i_test));

        end
    end
end


%% UL SINR
figure; 
Tau = 25;

y = reshape(UL_SINR_LS_ra,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'cx-');
hold on;
% 
y = reshape(UL_SINR_LS_sp,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'bx-');
hold on;


y = reshape(UL_SINR_LS_wg,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'gx-');
hold on;

y = reshape(UL_SINR_LS_sw,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'rx-');
hold on;

% y = reshape(UL_SINR_MMSE_ra,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'c*-');
% hold on;
% 
% y = reshape(UL_SINR_MMSE_sp,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'b*-');
% hold on;
% 
% 
% y = reshape(UL_SINR_MMSE_wg,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'g*-');
% hold on;
% 
% y = reshape(UL_SINR_MMSE_sw,1,L*K*Test_num);
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'r*-');
% hold on;

y = reshape(UL_SINR_MLE_ra,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'co-');
hold on;

y = reshape(UL_SINR_MLE_sp,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'bo-');
hold on;


y = reshape(UL_SINR_MLE_wg,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'go-');
hold on;

y = reshape(UL_SINR_MLE_sw,1,L*K*Test_num);
ymin=min(y); 
ymax=max(y); 
x=linspace(ymin,ymax,Tau); 
yy=hist(y,x);
yy=yy/length(y)/(x(2)-x(1));
s=0; 
for i=2:length(x) 
    s=[s,trapz(x([1:i]),yy([1:i]))]; 
end 
plot(x,s,'ro-');
hold on;


grid on;
%axis([-80 40 0 1]); 
title('CDF of uplink R');
% legend('LS:Random','LS:SPRS','LS:WGC-PA','LS:SPRS+WGC-PA','MMSE:Random','MMSE:SPRS','MMSE:WGC-PA','MMSE:SPRS+WGC-PA','MLE:Random','MLE:SPRS','MLE:WGC-PA','MLE:SPRS+WGC-PA');
legend('Random','SPRS','WGC-PA','SPRS+WGC-PA','MLE:Random','MLE:SPRS','MLE:WGC-PA','MLE:SPRS+WGC-PA');
xlabel('Average uplink achievable rate per user R (bps/Hz)');
ylabel('CDF');
hold off;

toc;















