% calcualte something via normal and IS and cmpare the marginal
% likelihoods.

N = 20;
beta = 2;
sigma = 1;
gamma = 1;
p = 0.8;

part = 1000;

Z0 = [2,2,2,0,1,0];  % initial state.
y = 4;         % number of observed events (over obs interval)
fs = 18;       % observed final size (over entire epidemic)

% particles in the initial state.
X0 = repmat(Z0,part,1);

%% Gillespie simulation


REPS1 = 1000;
mL_gillespie = zeros(REPS1,1);

for ii=1:REPS1

    % Forward simulation using Gillespie.
    X1 = SEEIIRp_gillespie_sim(N,beta,sigma,gamma,p,X0);
    
    % select particles that match.
    bum = X1(:,3) == Z0(3)+y &...
          X1(:,1) ~= X1(:,6) &...
          X1(:,1) <= fs &...
          X1(:,4) <= N-fs;

    mL_gillespie(ii) = sum(bum)/part;

end

like_gi = mean(mL_gillespie)


%% 

REPS2 = 100;
mL_is = zeros(REPS2,1);

for ii=1:REPS2
    % forwad sim using importance sampling.
    [X11,w] = SEEIIRp_is(N,beta,sigma,gamma,p,X0,y,fs);

    mL_is(ii) = mean(w);
end

like_is = mean(mL_is)
