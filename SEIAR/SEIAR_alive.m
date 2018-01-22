
function [X1,like] = SEIAR_alive(N,beta,sigma,gamma,p,X,NR,NF)
% given an initial state, simulate forward 1 day and return the state
% vector.

[part,~] = size(X); %number of particles
a = zeros(5,1);
X1 = zeros(part,5); % holds the particles that match.

jj = 1;
n = 0;    % counts number of times we need to simulate until we get part+1 matches
while jj <= part + 1
    
    % sample a particle
    ii = randi(part,1);
    n = n + 1;
    Z = X(ii,:); 
    
    y = Z(3) + NR;
    
    time = 0;

    while Z(1)-Z(5)-Z(4) > 0 && Z(3) <= y && Z(5) <= N-NF

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
    
    if Z(1)-Z(5)-Z(4) > 0 && Z(3) == y && Z(5) <= N-NF
        % particle matches observation
        % only add the first part matches and not the last.
        if jj <= part 
            X1(jj,:) = Z;
        end
        jj = jj+1;
    end
    
end

% has to have this form to be unbiased. 
like = part/(n-1);


