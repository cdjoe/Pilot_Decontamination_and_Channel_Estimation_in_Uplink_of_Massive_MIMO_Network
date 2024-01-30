function pilot = F_ZC(N, S)

% N = 10;   % the length of the pilot sequence (N>=K)
        pilot = zeros(N,S);   % the pilot signal
            
        if mod(N,2) == 0
            for k = 1:S
                for i = 0 : N-1
                    ll = mod(i+k-1, N);
                    pilot(i+1,k) =  exp(1i*pi*ll*ll/N);   % N is even.
                end
            end
        else
            for k = 1 : S
                for i = 0 : N-1
                    ll = mod(i+k-1, N);
                    pilot(i+1,k) = exp(1i*pi*ll*(ll+1)/N);   % N is odd.
                end
            end
        end
