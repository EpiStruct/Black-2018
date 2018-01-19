% This script estimates the likelihood for an SIR model using bootstrap
% sampling and importance sampling. 
% This is to demonstrate that it works correctly.

N = 20;
beta = 2;
gamma = 1;

y = 4;       % observed number of infections over a day.
X0 = [1,0];  % initial state.

%% Estimate the likelihood using the Gillespie algorithm

part = 1000;
REPS = 1000;

mL_gillespie = zeros(REPS,1);

for ii=1:REPS

    X1 = SIR_forward_day_gillespie(N,beta,gamma,repmat(X0,part,1));
    
    % particle mathces the data, y, and has not faded out.
    mL_gillespie(ii) = sum(X1(:,1)==X0(1)+y & X1(:,1)~=X1(:,2))/part;

end

like_gi = mean(mL_gillespie)


%% Estaimte the likelihood using importance sampling

part = 100; % use many less particles than above.
REPS = 100;

mL_is = zeros(REPS,1);

for ii=1:REPS
    % forwad sim using importance sampling.
    [X11,w] = SIR_forward_day_is(N,beta,gamma,repmat(X0,part,1),y);

    mL_is(ii) = mean(w);
end

like_is = mean(mL_is)

