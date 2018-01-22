% Gillespie simulation of the SEIR model.

function X = SEIR_gillespie_sim(N,beta,sigma,gamma,X)

[part,~] = size(X); % number of particles
a = zeros(3,1);

for kk = 1:part

    Z = X(kk,:); 
    time = 0;

    while Z(1)-Z(3) > 0 % E+I>0

        I = Z(2)-Z(3);
        E = Z(1) - Z(2);
        S = N - Z(1);
        
        a(1) = beta*S*I/(N-1);
        a(2) = sigma*E;
        a(3) = gamma*I;

        a0 = sum(a);

        % time to next event
        time = time - log(rand)/a0;

        if time > 1
            break;
        end

        r = rand;
        index = 0;
        tot = 0;
        while tot <= r*a0
            tot = tot +a(index+1);
            index = index +1;
        end

        Z(index) = Z(index) + 1;
    end
    
    X(kk,:) = Z;

end

end