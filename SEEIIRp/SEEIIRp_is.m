% SEEIIR model with partial observation of the E2->I1 event.

function [X,w] = SEEIIRp_is(N,beta,sigma,gamma,p,X,y,NF)

[part,~] = size(X); %number of rows
w = zeros(part,1);

a = zeros(6,1);
b = zeros(6,1);

% arrays to hold the stack.
t_next    = zeros(1,y+2);
type_next = zeros(1,y+2);

b = zeros(4,1);

for kk = 1:part

    Z = X(kk,:);
    
    type_next(1:y) = 3;
    type_next(y+1:y+2) = 8;
    
    r = y;
    
    % generate the order statistics. 
    ac = 0.0;
    for ii=1:y
        ac = ac - log(rand);
        t_next(y+1-ii) = ac;
    end  
    ac = ac - log(rand);
    
    for ii=y:-1:1
        t_next(ii) = t_next(ii)/ac;
    end

    L_imp = -log(factorial(y));
    t = 0;

    while r > 0
        
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
  
        % check if an event needs to forced within the next interval.
        if E2 == 0 && type_next(r) == 3
   
            if E1 == 0
                % force a type 1
                lambda = a(1);
                
                t_dash = trunc_exp(lambda,t_next(r)-t);
                wc = log(lambda) - lambda*t_dash - log(1-exp(-lambda*(t_next(r)-t)));
                L_imp = L_imp - wc;
                
                % add the time and event index to the stack.
                r = r + 1;
                type_next(r) = 1;
                t_next(r) = t + t_dash;
                
            else
                % force a type 2
                lambda = a(2);
                
                t_dash = trunc_exp(lambda,t_next(r)-t);
                wc = log(lambda) - lambda*t_dash - log(1-exp(-lambda*(t_next(r)-t)));
                L_imp = L_imp - wc;
                
                % add the time and event index to the stack.
                r = r + 1;
                type_next(r) = 2;
                t_next(r) = t + t_dash;
              
            end     
            
        end
        
        % set the rates of the modified process
        b(3) = 0;
        
        % if p==1 don't allow more Z1 than the final size.
        if type_next(r) == 1 || (Z(1)==NF && p ==1)
            b(1) = 0;
        else
            b(1) = a(1);
        end
        
        % if p==1 don't allow more Z1 than the final size.
        if type_next(r) == 2 
            b(2) = 0;
        else
            b(2) = a(2);
        end
        
        if Z(1)-Z(6) == 1 % E+I==1 to stop fadeout.
            b(6) = 0;
        else
            b(6) = a(6);
        end

        % stop too many non-detection events when we know the final size.
        if Z(4) == N-NF
            b(4) = 0;
        else
            b(4) = a(4);
        end
        
        b(5) = a(5);
        
        b0 = sum(a);
        a0 = sum(b);
        
        % propose a time for the next event.
        t_dash = -log(rand)/a0;
        
        
        if t_dash + t < t_next(r)
            % add an event at the proposed time.
            
            rn = rand;
            index = 0;
            tot =0;
            while tot <= rn*a0
                tot = tot +b(index+1);
                index = index +1;
            end
            
            orig = log(a(index))-b0*t_dash;
            is = log(b(index))-a0*t_dash;
            
            L_imp = L_imp + orig - is;
            
            t = t + t_dash;
            Z(index) = Z(index) +1;
            
            
        else
            
            ne = type_next(r); % which event are we adding?
            
            orig = log(a(ne))-b0*(t_next(r)-t);
            is = -a0*(t_next(r)-t);
            
            L_imp = L_imp + orig - is;

            t = t_next(r);
            Z(ne) = Z(ne) +1;
            
            % pop top element off the stack
            r = r -1; 
            
        end

    end

    % simulate to the end of the day. 

    while 1==1

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
               
        % set the rates of the modified process
        b = a;
        
        b(3) = 0;
        
        % if p==1 don't allow more Z1 than the final size.
        if Z(1)==NF && p ==1
            b(1) = 0;
        else
            b(1) = a(1);
        end
        
    
        if Z(1)-Z(6) == 1 % E+I==1 to stop fadeout.
            b(6) = 0;
        else
            b(6) = a(6);
        end

        % stop too many non-detection events when we know the final size.
        if Z(4) == N-NF
            b(4) = 0;
        else
            b(4) = a(4);
        end

        b0 = sum(a);
        a0 = sum(b);
        
        % propose a time for the next event.
        t_dash = -log(rand)/a0;
        
        
        if t_dash + t < 1
            % add an event at the proposed time.
            
            rn = rand;
            index = 0;
            tot =0;
            while tot <= rn*a0
                tot = tot +b(index+1);
                index = index +1;
            end
            
            orig = log(a(index))-b0*t_dash;
            is = log(b(index))-a0*t_dash;
            
            L_imp = L_imp + orig - is;
            
            t = t + t_dash;
            Z(index) = Z(index) +1;
            
            
        else
            
            % end of the day reached.
            orig = -b0*(1-t);
            is = -a0*(1-t);

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


