%% Symmetrical NMF

function H = SNMF(G, k, MaxIter, epsilon)

%%% Input
%
% G: the adjacency matrix of a network
% k: the number of clusters
% MaxIter: the maximal number of iterations for alternating minimization
% epsilon: the convergence parameter

%% Normalize the network

G = G/sqrt(trace(G'*G));

%% Initializaiton

g = size(G,1);
H = rand(g,k);
H = H/sqrt(trace(H'*H));

J1 = norm(G - H*H', 'fro')^2;
Iter = 1;
beta = 0.5;
Delta = 99999;

%% Alternating update

while Delta > epsilon && Iter <= MaxIter
    
    Num = G*H;
    Denom = (H')*H;
    Denom = H*Denom;
    H = H.*(1 - beta + beta*(Num./Denom));
    
    J2 = norm(G - H*H', 'fro')^2;
    Delta = J1 - J2;
    J1 = J2;
    Iter = Iter + 1;
    
end

end