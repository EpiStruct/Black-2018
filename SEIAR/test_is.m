% Calcualting the likelihood using Gillespie simulaions, the alive filter
% and importance sampling.

% Note that the number of particles and repetitions are arbitrarily set
% for each algorithm.


N = 164;
beta = 2;
sigma = 1;
gamma = 1;
p = 0.9;

Z0 = [1,1,0,0,0];  % initial state.
y = 4;             % number of observed events (over obs interval)
fs = 146;          % observed final size (over entire epidemic)

%% simulate using basic Gillespie algorithm

part = 1000;
X0 = repmat(Z0,part,1);
REPS1 = 500;

mL_gillespie = zeros(REPS1,1);

for ii=1:REPS1
  
    X1 = SEIAR_gillespie_sim(N,beta,sigma,gamma,p,X0);
    
    % find particles that match.
    match = X1(:,3) == Z0(3)+y &...
            X1(:,1) ~= X1(:,5) + X1(:,4) &...
            X1(:,5) <= N-fs;

    mL_gillespie(ii) = sum(match)/part;
   
end

like_gi = mean(mL_gillespie)



%% Simulate using the alive filter

part = 100;
X0 = repmat(Z0,part,1);

REPS1 = 10;

mL_alive = zeros(REPS1,1);

for ii=1:REPS1
  
    [X1,like] = SEIAR_alive(N,beta,sigma,gamma,p,X0,y,fs);
    mL_alive(ii) = like;
    
end

like_al = mean(mL_alive)

%% Simulate using importance sampling

REPS2 = 500;
part = 100;
mL_is = zeros(REPS2,1);
X0 = repmat(Z0,part,1);

for ii=1:REPS2

    [X11,w] = SEIAR_is(N,beta,sigma,gamma,p,X0,y,fs);
    mL_is(ii) = mean(w);

end

like_is = mean(mL_is)

