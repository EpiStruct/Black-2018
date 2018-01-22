% Gillespie simulation of the SEEIIR model with partial detection.

function X = SEEIIRp_gillespie_sim(N,beta,sigma,gamma,p,X)

[part,~] = size(X); % number of particles
a = zeros(6,1);

for kk = 1:part

    Z = X(kk,:); 
    time = 0;

    while Z(1)-Z(6) > 0 % E1+E2+I1+I2 > 0

        S = N - Z(1);
        E1 = Z(1) - Z(2);
        E2 = Z(2) - Z(3) - Z(4);
        I1 = Z(3) + Z(4) - Z(5);
        I2 = Z(5) - Z(6);
        
        a(1) = beta*S*(I1+I1)/(N-1);
        a(2) = 2*sigma*E1;
        a(3) = 2*sigma*E2*p;
        a(4) = 2*sigma*E2*(1-p);
        a(5) = 2*gamma*I1;
        a(6) = 2*gamma*I2;
        
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