clear;

clc; 
tic;
%% system parameters
M = 512;   
L = 19;
K = 10;
S = 15;

R_Cell = 1000;
r_Min = 100;
alpha = 3.8;
sigma_shadow = 8;
rho_ul = 10^(15/10);                                                        % 15 dB                                                     % 15 dB
rho_pilot = 10^(15/10);

g = 20;   % number of grids
Lambda = zeros(1,g);

% First iteration
Max = 1;   
min = 0.05;
delta = (Max-min)/(g-1);

% Second iteration
% Max = 0.4 + (0.05/2);      % lambda1_max+(delta1/2)
% min = 0.4 - (0.05/2);            % lambda1_max-(delta1/2)
% delta = (Max-min)/(g-1);

Lambda(1,1) = min;
for i = 1:g-1
    Lambda(1,i+1) = Lambda(1,i)+delta;
end

% Test_num = 300;
Test_num = 50;
UL_SINR_CS_MF = zeros(K, L, Test_num);
UL_SINR_LS_sp = zeros(K, L, Test_num, g);
UL_SINR_WGCPA_MF = zeros(K, L, Test_num);
UL_SINR_WGCPAsoft_MF = zeros(K, L, Test_num);
SS = zeros(Test_num,g);

%% Simulation
for i_test = 1:Test_num
    display(i_test);
    
    for j = 1:g
        lambda = Lambda(1,j);
%         display(lambda);


            %% Generate Channel vector
            [H, Beta] = F_H_Generate(M, L, K, R_Cell, r_Min, sigma_shadow, alpha);

            %% Random Pilot Assignment: k-th pilot for k-th user
    %         P = zeros(L, K);
    %         for i = 1:L                                                         % random pilot assignment
    %             temp = randperm(S);
    %             P(i,:) = temp(1:K);
    %         end
    %         
    % %         R_p = P.' * conj(P);
    %         
    %         UL_SINR_CS_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);



            %% Random Pilot Assignment For Certen_SOFT         
%             [pilots, N] = F_SPRS2(L,K,Beta,lambda);
% 
%             Y = received_signal(M,N,L,rho_ul,H,pilots);
% 
%             % ===== LS --> SINR =====%
%             H_est = CE_LS(M, K, L, Y, pilots, rho_ul, N);  
% 
%             UL_SINR_LS_sp(:, :, i_test, j) = F_SINR(L,K,H,H_est,rho_ul,pilots);

    %         UL_SINR_CSsoft_MF(:,:,i_test,j) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);
            [P, eu]= F_CS_SOFT(L, K, S, Beta,lambda);

        
            UL_SINR_LS_sp(:,:,i_S,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1); 

            %% Weighted Graph Coloring Based Pilot Assignment WGC-PA
    %         P = F_WGCPA_Pilot(L, K, S, Beta);
    %         
    %         UL_SINR_WGCPA_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P, rho_pilot, rho_ul, 1);


             %% Weighted Graph Coloring Based Pilot Assignment For Certen_SOFT
    %         [P1 eu]= F_WGCPA_soft_Pilot(L, K, S, Beta,lambda);
    %         
    %         UL_SINR_WGCPAsoft_MF(:,:,i_test) = F_ULSINR_CAL(L, K, M, H, P1, rho_pilot, rho_ul, 1);

            SS(i_test,j) = SS(i_test,j) + N;
    end

end

%% UL SINR
figure; 

% y = 10*log10(reshape(UL_SINR_CS_MF,1,L*K*Test_num));
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'co-');
% hold on;
% 
% y = 10*log10(reshape(UL_SINR_CSsoft_MF,1,L*K*Test_num));
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'bo-');
% hold on;
% 
% 
% y = 10*log10(reshape(UL_SINR_WGCPA_MF,1,L*K*Test_num));
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'go-');
% hold on;
% 
% y = 10*log10(reshape(UL_SINR_WGCPAsoft_MF,1,L*K*Test_num));
% ymin=min(y); 
% ymax=max(y); 
% x=linspace(ymin,ymax,Tau); 
% yy=hist(y,x);
% yy=yy/length(y)/(x(2)-x(1));
% s=0; 
% for i=2:length(x) 
%     s=[s,trapz(x([1:i]),yy([1:i]))]; 
% end 
% plot(x,s,'ro-');
% hold on;

Bw = 1;
mu0=0.05;
T_1=1-S/K*mu0;


SINR = zeros(g, 1);
R = zeros(g, 1);
for j = 1:g
    for i_test = 1:Test_num
        T_soft=1-SS(i_test,j)/K*mu0;
        
        for i = 1:L
            for k = 1:K
                SINR(j) = SINR(j) + UL_SINR_LS_sp(k,i,i_test,j)/Test_num/L/K;
%                 R(j) = R(j) + T_soft*log2(1+UL_SINR_LS_sp(k,i,i_test,j))/Test_num/L/K;
            end
        end
    end
end
% plot(Lambda,R,'g-o');
% grid on;
% % axis([0.375 0.425 7.34 7.4]);
% title('The IGS procedure to obtain the near-optimal adjustment parameter \lambda');
% % legend('Random','SPRS','WGC-PA','SPRS+WGC-PA', 'Location', 'Best');
% legend('First iteration of the IGS');
% % legend('Second iteration of the IGS');   
% xlabel('\lambda');
% ylabel('Average uplink achievable rate per user R (bps/Hz)');
% hold off;

%  
plot(Lambda,SINR,'r-o');
grid on;
%axis([-80 40 0 1]);
title('The IGS procedure to obtain the near-optimal adjustment parameter \lambda');
% legend('Random','SPRS','WGC-PA','SPRS+WGC-PA', 'Location', 'Best');
legend('First iteration of the IGS');
% % legend('Second iteration of the IGS');   
xlabel('\lambda');
ylabel('SINR(dB)');
hold off;






toc;















