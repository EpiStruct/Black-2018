% SEEIIR model with partial observation of E->I event.

function [X,w] = SEIAR_is(N,beta,sigma,gamma,p,X,y,NF)

[part,~] = size(X); %number of rows
w = zeros(part,1);

% arrays to hold the stack.
t_next    = zeros(1,y+2);
type_next = zeros(1,y+2);

b = zeros(5,1);
a = zeros(5,1);

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

    L_imp = -gammaln(y+1);
    t = 0;

    while r > 0
        
        S = N - Z(1);
        E = Z(1) - Z(2) - Z(5);
        Ip = Z(2) - Z(3);
        Is = Z(3) - Z(4);
                
        a(1) = beta*S*(Ip+Is)/(N-1);
        a(2) = sigma*p*E;
        a(3) = gamma*Ip;
        a(4) = gamma*Is;
        a(5) = sigma*(1-p)*E;

        if type_next(r) == 3 && Ip == 0
            % Next event is 3, but Ip=0
            
            if E == 0
                % E = 0, so need to put int an type 1 first.
                
                lambda = a(1);
                t1 = log(1-rand*(1-exp(-lambda*(t_next(r)-t))))/(-lambda);

                wc = log(lambda) - lambda*t1 - log(1-exp(-lambda*(t_next(r)-t)));
                L_imp = L_imp - wc;
                
                r = r + 1;
                type_next(r) = 1;
                t_next(r) = t + t1;
                
            else
                % E > 0, so can put in the type 2 now.
                
                lambda = a(2);
                t1 = log(1-rand*(1-exp(-lambda*(t_next(r)-t))))/(-lambda);
                wc = log(lambda) - lambda*t1 - log(1-exp(-lambda*(t_next(r)-t)));
                L_imp = L_imp - wc;
                
                r = r + 1;
                type_next(r) = 2;
                t_next(r) = t + t1;
                
            end
            
            
        elseif type_next(r) == 2 && E == 0
            % next event is a type 2, but E = 0, so force a type 1.
            
            lambda = a(1);
            t1 = log(1-rand*(1-exp(-lambda*(t_next(r)-t))))/(-lambda);
            wc = log(lambda) - lambda*t1 - log(1-exp(-lambda*(t_next(r)-t)));
            L_imp = L_imp - wc;
            
            r = r + 1;
            type_next(r) = 1;
            t_next(r) = t + t1;
            
        end
        
        
        % set the rates of the modified process
        b(3) = 0;
        
        % if p==1 don't allow more Z1 than the final size.
        if type_next(r) == 1 || (Z(1)==NF && p ==1)
            b(1) = 0;
        else
            b(1) = a(1);
        end
        
        % no more Z2 if we have reached the final size.
        if type_next(r) == 2 ||  Z(2)==NF
            b(2) = 0;
        else
            b(2) = a(2);
        end
        
   
        
        if Z(1)-Z(4)-Z(5) == 1 % E+Ip+Is == 1 to stop fadeout.
            b(4) = 0;
            b(5) = 0;
        else
            b(4) = a(4);
           
            % stop too many non-detection events when we know the final size.
            if Z(5) == N-NF
                b(5) = 0;
            else
                b(5) = a(5);
            end

        end

        a0 = a(1)+a(2)+a(3)+a(4)+a(5);
        b0 = b(1)+b(2)+b(4)+b(5);
        
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

    % after all the detection events have been put in we have to simulate to
    % the end of the day. 

    while 1==1

        S = N - Z(1);
        E = Z(1) - Z(2) - Z(5);
        Ip = Z(2) - Z(3);
        Is = Z(3) - Z(4);
                
        a(1) = beta*S*(Ip+Is)/(N-1);
        a(2) = sigma*p*E;
        a(3) = gamma*Ip;
        a(4) = gamma*Is;
        a(5) = sigma*(1-p)*E; 
         

        % set the rates of the modified process
        b(3) = 0;
        
        % if p==1 don't allow more Z1 than the final size.
        if (Z(1)==NF && p ==1)
            b(1) = 0;
        else
            b(1) = a(1);
        end
        
        % no more Z2 if we have reached the final size.
        if Z(2)==NF
            b(2) = 0;
        else
            b(2) = a(2);
        end

        
        if Z(1)-Z(4)-Z(5) == 1 % E+Ip+Is == 1 to stop fadeout.
            b(4) = 0;
            b(5) = 0;
        else
            b(4) = a(4);
           
            % stop too many non-detection events when we know the final size.
            if Z(5) == N-NF
                b(5) = 0;
            else
                b(5) = a(5);
            end

        end

        a0 = a(1)+a(2)+a(3)+a(4)+a(5);
        b0 = b(1)+b(2)+b(4)+b(5);
        
        % propose a time for the next event.
        t_dash = -log(rand)/b0;
        
        
        if t_dash + t < 1
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
            
            % end of the day.
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


