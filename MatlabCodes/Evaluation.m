%% Evaluate clustering accuracy and NMI

function [Accs, AvgAcc, AvgNMI] = Evaluation(Us, DomLabels)

%%% Input
%
% Us: a set of cluster indicator matrices
% DomLabels: the ground truth of domain clusters

%% Initialization

g = length(Us);
Corrects = zeros(g,1);
Totals = zeros(g,1);
NMIs = zeros(g,1);

%% Accuracy and NMI

for i = 1:g
    
    U_i = Us{i};
    t_i = size(U_i,2);
    DomLabel_i = DomLabels{i};
    [Vals, U_i_idx] = max(U_i,[],2);
    
    NonNoise_idx = DomLabel_i > 0; % Eliminate noisy nodes that have label 0
    DomLabel_i = DomLabel_i(NonNoise_idx);
    U_i_idx = U_i_idx(NonNoise_idx);
    
    for j = 1:t_i
        
        ClusLabel_j = DomLabel_i(U_i_idx == j);
        [Maj, Freq] = mode(ClusLabel_j);
        Corrects(i) = Corrects(i) + Freq;
        Totals(i) = Totals(i) + length(ClusLabel_j);
        
    end
    
    NMI_u_i = NMI(DomLabel_i, U_i_idx);
    NMIs(i) = NMI_u_i;
    
end

Accs = Corrects./Totals; % Accuracies of individual domain-specific networks
AvgAcc = sum(Corrects)/sum(Totals); % The overall accuracy
AvgNMI = mean(NMIs); % The average NMI

end