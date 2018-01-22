% importance sampling for the SEIR model with observation of E->I event.

function [X,w] = SEIR_one_day_is(N,beta,sigma,gamma,X,NR,NF)

[part,~] = size(X);     % number of particles
w = zeros(part,1);

% arrays to hold the stack.
t_next    = zeros(1,NR+1);
type_next = zeros(1,NR+1);

b = zeros(3,1);
a = zeros(3,1);

for kk = 1:part

    Z = X(kk,:);
    
    type_next(1:NR) = 2;
    type_next(NR+1) = 8;
    
    r = NR;
    
    % generate the order statistics. 
    ac = 0.0;
    for ii=1:NR
        ac = ac - log(rand);
        t_next(NR+1-ii) = ac;
    end  
    ac = ac - log(rand);
    
    for ii=NR:-1:1
        t_next(ii) = t_next(ii)/ac;
    end

    L_imp = -log(factorial(NR));
    t = 0;

    while r > 0
        

        I = Z(2)-Z(3);
        E = Z(1) - Z(2);
        S = N - Z(1);
        
        % calculate the rates of the original process.
        a(1) = beta*S*I/(N-1);
        a(2) = sigma*E;
        a(3) = gamma*I;
  
        % check to see if e=1 events need to to be forced within the next
        % interval.
        if E == 0 && type_next(r)==2
                            
            % generate a time given the current state
            % using a truncated exponential RV.
            lambda = a(1);
            
            t_dash = trunc_exp(lambda,t_next(r)-t);
            % is contribution from this is...
            wc = log(lambda) - lambda*t_dash - log(1-exp(-lambda*(t_next(r)-t)));
            L_imp = L_imp - wc;
            
            % add the time and event index to the stack.
            r = r + 1;
            type_next(r) = 1;
            t_next(r) = t + t_dash;
            
        end
        
        % set the rates of the modified process
        b(2) = 0;
        
        if type_next(r) == 1 || Z(1)==NF
            b(1) = 0;
        else
            b(1) = a(1);
        end
        
        if Z(1)-Z(3) == 1 % E+I==1 to stop fadeout.
            b(3) = 0;
        else
            b(3) = a(3);
        end

        a0 = sum(a);
        b0 = sum(b);
        
        % propose a time for the next event.
        t_dash = -log(rand)/b0;
        
        
        if t_dash + t < t_next(r)
            % add an event at the proposed time.
            
            rn = rand;
            index = 0;
            tot =0;
            while tot <= rn*b0
                tot = tot +b(index+1);
                index = index +1;
            end
            
            orig = log(a(index))-a0*t_dash;
            is = log(b(index))-b0*t_dash;
            
            L_imp = L_imp + orig - is;
            
            t = t + t_dash;
            Z(index) = Z(index) +1;
            
            
        else
            
            ne = type_next(r); % which event are we adding?
            
            orig = log(a(ne))-a0*(t_next(r)-t);
            is = -b0*(t_next(r)-t);
            
            L_imp = L_imp + orig - is;

            t = t_next(r);
            Z(ne) = Z(ne) +1;
            
            % pop top element off the stack
            r = r -1; 
            
        end

    end

    % after all the detection events (en=2) have been put in we still need 
    % to simulate to the end of the day. 

    while 1==1

        I = Z(2)-Z(3);
        E = Z(1) - Z(2);
        S = N - Z(1);
        
        a(1) = beta*S*I/(N-1);
        a(2) = sigma*E;
        a(3) = gamma*I;

        a0 = sum(a);
        
        b = a;
        b(2) = 0;
        
        % limit the number of Z1 events.
        if Z(1)==NF
            b(1) = 0;
        end
        
        if Z(1)-Z(3) == 1
            b(3) = 0;
        end
        
        b0 = sum(b);

        t_dash = -log(rand)/b0;

        if t_dash + t < 1
            % then we add

            rn = rand;
            index = 0;
            tot =0;
            while tot <= rn*b0
                tot = tot +b(index+1);
                index = index +1;
            end

            orig = log(a(index))-a0*t_dash;
            is = log(b(index))-b0*t_dash; 

            L_imp = L_imp + orig - is;

            t = t + t_dash;
            Z(index) = Z(index) +1;

        else
            % move to the end of the day.
            orig = -a0*(1-t);
            is = -b0*(1-t);

            L_imp = L_imp + orig - is;

            break;

        end
    
    end

    w(kk) = exp(L_imp);
    X(kk,:) = Z;
    
end    
    
end

function po = trunc_exp(lambda,t) 

    po = log(1-rand(1)*(1-exp(-lambda*t)))/(-lambda);

end
