% Code for the decay model example in section 2.

part = 10^5;
gamma = 1;
y = 10;


%% Simulations using basic Gillespie algorithm.
X = zeros(part,1);

for ii=1:part
    
    Z = X(ii);
    t = 0;
    
    while 1==1
    
        % rate of event.
        a = gamma*(20-Z);

        % time until next event is
        t_dash = -log(rand)/a;

        t = t + t_dash;

        if t > 1
            break;

        end

        Z = Z + 1;
    
    end
    
    X(ii) = Z;
    
end

% likelihood is estimated from the number of matches.
like_gl = sum(X==y)/part


%% now do the same using importance sampling

X = zeros(part,1);
w = zeros(part,1);

NR = 10;

for ii=1:part
    
    Z = X(ii);
    t = 0;
    
    % generate the times.
    t_next = sort(rand(NR,1));
    n = 1; 
    
    % contribution to the log importnace weight from the forced events.
    L_imp = -log(factorial(NR));
    
    while n <= NR
    
        % the rate of some event.
        a = gamma*(20-Z);

        % transition density of next event at time t_next(n)

        L_imp = L_imp - a*(t_next(n)-t) + log(a);
        
        t = t_next(n);
        n = n + 1;
        Z = Z + 1;

    end
    
    % now calcualte the final contribution, the probability of no more
    % events by the end of the day.
    a = gamma*(20-Z);

    L_imp = L_imp - a*(1-t);

    X(ii) = Z;
    w(ii) = L_imp;
    
end

like_is = mean(exp(w))

