% simulates forward 1 day using Gillespie algorithm

function X = forward_day_gillespie(N,beta,gamma,X)

    [part,~] = size(X); %number of particles
    b = zeros(2,1);
    beta_norm = beta/(N-1);
    
    for kk = 1:part

        Z = X(kk,:);  
        t = 0;

        while Z(1)-Z(2) > 0

            b(1) = beta_norm*(N-Z(1))*(Z(1)-Z(2));
            b(2) = gamma*(Z(1)-Z(2));
            
            b0 = b(1)+b(2);

            t_dash = -log(rand)/b0;

            if t + t_dash > 1
                % gone over the next day
                break;
            end

            t = t+t_dash;

            if rand*b0<b(1)
                Z(1) = Z(1) +1;
            else
                Z(2) = Z(2) +1;
            end

        end
        
        X(kk,:) = Z;
    end
end


