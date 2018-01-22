% Gillespie simulation of the SEIAR model.

function X = SEIAR_gillespie_sim(N,beta,sigma,gamma,p,X)

[part,~] = size(X); % number of particles
a = zeros(5,1);

for kk = 1:part

    Z = X(kk,:); 
    time = 0;

    while Z(1)-Z(5)-Z(4) > 0 % E+Ip+Is>0

        S = N - Z(1);
        E = Z(1) - Z(2) - Z(5);
        Ip = Z(2) - Z(3);
        Is = Z(3) - Z(4);
                
        a(1) = beta*S*(Ip+Is)/(N-1);
        a(2) = sigma*p*E;
        a(3) = gamma*Ip;
        a(4) = gamma*Is;
        a(5) = sigma*(1-p)*E;
        
        a0 = a(1)+a(2)+a(3)+a(4)+a(5);

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