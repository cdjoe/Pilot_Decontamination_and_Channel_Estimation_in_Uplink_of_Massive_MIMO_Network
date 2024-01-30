function [pilots, N] = F_SPRS_WGCPA2(L, K, Beta, lambda, S)

%===== SPRS =====%
phi = zeros(1,L);
cent = zeros(K,L); 
edge = zeros(K,L); 
beta_2=zeros(K,L);
sum_c = zeros(1,L);
% phi = lambda/K;

for l=1:L
    for k=1:K
%    sum = sum + Beta(k,l,l)^2;
        beta_2(k,l) = Beta(k,l,l)^2;
    end 
    phi(1,l)=sum(beta_2(:,l))*lambda/K;
end

for l = 1:L
    for k = 1:K
        if phi(1,l) < beta_2(k,l)     
           cent(k,l) = cent(k,l) + 1;   % the number of center users
           sum_c(1,l) = sum_c(1,l) + 1;
        else
           edge(k,l) = edge(k,l) + 1;       % the number of edge users
        end
    end
end
            
c = S;     % max(sum_c)
e = sum(edge,'all');
SS = c + e;
N = SS;

%===== ZC =====%
pilot = F_ZC(N,SS);
pilots_c = pilot(:,e+1:SS);
pilots_e = pilot(:,1:e);

pilots = zeros(N,K,L);
n = 1;
for i = 1:L   
    for k = 1:K       
        if edge(k,i) == 1
            pilots(:,k,i)=pilots_e(:,n);
            n = n + 1;
        end
    end
end

%===== WGC-PD =====%
eta = zeros(L,K,L,K);
for i1 = 1:L                                                            
    for k1 = 1:K        % pilot contamination strength calculation
        for i2 = 1:L
            for k2 = 1:K
                if (i1~=i2)
                    eta(i1,k1,i2,k2) = Beta(k1,i2,i1)^2/Beta(k2,i2,i2)^2+Beta(k2,i1,i2)^2/Beta(k1,i1,i1)^2;
                end
            end
        end
    end
end

eta_max = 0;
for j_bs = 1:L                    % First 2 users pilot assignment                                 
    for k_user = 1:K
        if pilots(1,k_user,j_bs) == 0
            for j1_bs = 1:L
                for k1_user = 1:K
                    if pilots(1,k1_user,j1_bs) == 0
                        if j_bs~=j1_bs && eta(j_bs,k_user,j1_bs,k1_user)>eta_max
                            j_2 = j_bs;
                            k_2 = k_user;
                            j_1 = j1_bs;
                            k_1 = k1_user;
                            eta_max = eta(j_bs,k_user,j1_bs,k1_user);
                        end
                    end
                end
            end
        end
    end
end
pilots(:,k_1,j_1) = pilots_c(:,1);
pilots(:,k_2,j_2) = pilots_c(:,2);

for t = (3+e):L*K                   % Pilot assignment for other users                                                        
    delta = zeros(L,K);                             % UE selection                     
    delta_max = 0;
    for j = 1:L
        for k = 1:K
            if pilots(1,k,j) == 0
               for j1 = 1:L
                    for k1 = 1:K
                        if j~=j1 && pilots(1,k1,j1) ~= 0
                            delta(j,k) = delta(j,k) + eta(j,k,j1,k1);
                        end
                    end
               end
               if delta(j,k)>delta_max
                    j_0 = j;
                    k_0 = k;
                    delta_max = delta(j,k);
               end
            end
        end
    end

    C = zeros(1,c);                              % Pilot selection in j0-cell 
    for k = 1:K
%                 if pilots(1,k,j_0) == 0
            for s = 1:c
                if pilots(:,k,j_0) == pilots_c(:,s)
                    C(s) = 1;
                end
            end
%                 end
    end

    zeta = zeros(1,c);
    for j = 1:L
        for k = 1:K
            for s = 1:c
                if pilots(:,k,j) == pilots_c(:,s)
                    zeta(s) = zeta(s) + eta(j,k,j_0,k_0);
                end
            end
        end
    end

    zeta_min = 1000;
    c0 = 0;
    for k = 1:c
        if zeta(k)<zeta_min && C(k)==0
            zeta_min = zeta(k);
            c0 = k;
        end
    end           
    pilots(:,k_0,j_0) = pilots_c(:,c0);

end

end