function [pilots, N] = F_SPRS2(L, K, Beta, lambda, S)

% L:cell 
% K:user 
% Beta:LSF coefficient 
% lambda:adjustment parameter
% N:length
% S:number of pilot

phi = zeros(1,L);   % threshold
cent = zeros(K,L); sum_c = zeros(1,L);
edge = zeros(K,L); 
beta_2=zeros(K,L);
%         phi = lambda/K;

for l=1:L
    for k=1:K
%                sum = sum + Beta(k,l,l)^2;
        beta_2(k,l) = Beta(k,l,l)^2;
    end

    phi(l)=sum(beta_2(:,l))*lambda/K;

end

for l = 1:L
    for k = 1:K
        if phi < beta_2(k,l)     
           cent(k,l) = cent(k,l) + 1;     % the center users
           sum_c(1,l) = sum_c(1,l) + 1;
        else
           edge(k,l) = edge(k,l) + 1;     % the edge users
        end
    end
end

c = max(sum_c);     % max(sum_c)
e = sum(edge,'all');
SS = c + e;

N = SS;

%=====ZC=====%
pilot = F_ZC(N,SS);
pilot_c = pilot(:,e+1:SS);
pilot_e = pilot(:,1:e);
temp = zeros(L,c);
for i=1:L
    temp(i,:) = randperm(c);
end

%===== Pilot Assignment =====%
n = 1; 
pilots = zeros(N,K,L);
for i = 1:L
    m = 1;
    for k = 1:K       
        if edge(k,i) == 1
            pilots(:,k,i)=pilot_e(:,n);
            n = n + 1;
        else
            pilots(:,k,i)=pilot_c(:,temp(i,m));   %m
            m = m + 1;
        end

    end
end
        
end

