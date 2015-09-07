%% Run NoNClus on the 20Newsgroup dataset

function Run20Newsgroup(a, k, t_u, t_v, MaxIter, epsilon, ComRate, Runtime)

%%% Input
%
% If no input parameters are provided, the default values will be used.
%
% a: a regularization parameter of NoNClus
% k: the number of main clusters
% t_u: a vector of the numbers of domain clusters in domain-specific networks
% t_v: a vector of the numbers of domain clusters in hidden factor matrices
% MaxIter: the maximal number of iterations for alternating minimization
% epsilon: the convergence parameter
% ComRate: the common node ratio between domain-specific networks
% Runtime: the run time on the 20Newsgroup dataset

%% Parameter initialization

g = 30; % The number of domain-specific networks

if ~exist('a', 'var') || isempty(a)
    a = 1;
end
if ~exist('k', 'var') || isempty(k)
    k = 3;
end
if ~exist('t_u', 'var') || isempty(t_u)
    t_u = 4*ones(g,1);
end
if ~exist('t_v', 'var') || isempty(t_v)
    t_v = 4*ones(k,1);
end
if ~exist('MaxIter', 'var') || isempty(MaxIter)
    MaxIter = 1000;
end
if ~exist('epsilon', 'var') || isempty(epsilon)
    epsilon = 1e-6;
end
if ~exist('ComRate', 'var') || isempty(ComRate)
    ComRate = 0.8;
end
if ~exist('Runtime', 'var') || isempty(Runtime)
    Runtime = 1;
end

AllAvgAcc = zeros(1,Runtime);
AllAccs = zeros(g,Runtime);
AllAvgNMI = zeros(1,Runtime);

%% Computation

for i = 1:Runtime
    
    % Generate NoN
    
    [DomNets, DomIDs, DomLabels, MainNet] = GenNoN_20Newsgroup(ComRate);
    
    % NoNClus
    
    Us = NoNClus(DomNets, DomIDs, MainNet, a, k, t_u, t_v, MaxIter, epsilon);
    
    % Evaluation
    
    [Accs, AvgAcc, AvgNMI] = Evaluation(Us, DomLabels);
    
    AllAvgAcc(i) = AvgAcc;
    AllAccs(:,i) = Accs;
    AllAvgNMI(i) = AvgNMI;
    
end

%% Results

disp('=============Results: NoNClus===============')

disp(['Overall accuracy: ' num2str(mean(AllAvgAcc)) ' +- ' num2str(std(AllAvgAcc))])
disp(['Average NMI: ' num2str(mean(AllAvgNMI)) ' +- ' num2str(std(AllAvgNMI))])

disp('============================================')

end