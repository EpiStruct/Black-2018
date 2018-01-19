% Importance sampling for the SIR model.

% Note that this code has a lot of redundnacy in it, for example introducing
% the variables a and a0. This is to make the comparison with the more
% complex models easier. 

function [X,omega] = forward_day_is(N,beta,gamma,X,NR)


    [part,~] = size(X); %number of rows
    omega = zeros(part,1);
    b = zeros(2,1);
    a = zeros(2,1);
    
    for kk = 1:part

        Z = X(kk,:);    

        % generate the exact infection times
        tR = sort(rand(1,NR));

        r = 1;
        L_imp = -gammaln(NR+1);
        t = 0;
        
        while r <= NR
            
            S = N-Z(1);
            I = Z(1)-Z(2);
            
            b(1) = beta*S*I/(N-1);
            b(2) = gamma*I;

            if I > 1
                a(2) = b(2);
            else
                a(2) = 0;
            end
            
            b0 = b(1)+b(2);
            a0 = a(1)+a(2);
            
            t_dash = -log(rand)/a0;

               if t_dash + t < tR(r)
                   % put the event in.

                    Z(2) = Z(2) +1;
                    t = t + t_dash;

                    orig = log(b(2))-b0*(t_dash);

                    % prob density of recovery event at t_dash.
                    is = log(a0) - a0*t_dash;

                    L_imp = L_imp + orig - is;
               else
                   % next event is an infection
                   orig = log(b(1))-b0*(tR(r)-t);

                   % is is the prob of not putting in an infection or that
                   % t_dash > tR(r)-t.

                   is = -a0*(tR(r)-t); % log prob of no event in the interval.

                   L_imp = L_imp + orig - is;

                   Z(1) = Z(1) +1;
                   t = tR(r);     
                   r = r+1;

               end

        end

        % after the end of all the infection events then er need to put in the
        % last recoveries.
        % in this part we can't let the disease fade out.

        while t < 1
                        
            S = N-Z(1);
            I = Z(1)-Z(2);
            
            b(1) = beta*S*I/(N-1);
            b(2) = gamma*I;
             
            if I > 1
                a(2) = b(2);
            else
                a(2) = 0;
            end
            
            b0 = b(1)+b(2);
            a0 = a(1)+a(2);
            
            t_dash = -log(rand)/a0;
            
            if t + t_dash > 1
                % gone over the next day
                
                orig = -b0*(1-t);
                is = -a0*(1-t);
                
                L_imp = L_imp + orig - is;
                break;
                
            else
                % still within the day
                Z(2) = Z(2) +1;
                t = t + t_dash;
                
                orig = log(b(2)) -b0*t_dash;
                is = log(a(2)) - a0*t_dash;
                
                L_imp = L_imp + orig - is;
            end

        end

        omega(kk) = exp(L_imp);
        X(kk,:) = Z;
        
    end

end


