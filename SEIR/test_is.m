% Test of importance sampling against basic gillespie simulations.

addpath ..

N = 20;
beta = 2;
sigma = 2;
gamma = 1;

part = 500;

Z0 = [1,0,0];  % initial state. 
y = 5;         % number of observed events (over interval [0,1] )
fs = 18;       % observed final size (over entire epidemic)

fprintf('Initial state: S=%d, E=%d, I=%d\n',N-Z0(1),Z0(1)-Z0(2),Z0(2)-Z0(3));

X0 = repmat(Z0,part,1);

%% Forward simulation using Gillespie algorithm.

reps = 500;
mL_gillespie = zeros(reps,1);

for ii=1:reps

    X1 = SEIR_gillespie_sim(N,beta,sigma,gamma,X0);
    mL_gillespie(ii) = sum(X1(:,2)==Z0(2)+y & X1(:,1) ~= X1(:,3) & X1(:,1) <=fs )/part;
    
end

like_gi = mean(mL_gillespie)


%% forward sim using importance sampling.

reps = 100;
mL_is = zeros(reps,1);

tic
for ii=1:reps

    [X11,w] = SEIR_is(N,beta,sigma,gamma,X0,y,fs);
    mL_is(ii) = mean(w);
end
toc

like_is = mean(mL_is)


