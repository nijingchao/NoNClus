%% Generate domain-specific networks from the 20Newsgroup dataset

function [DomNets, DomIDs, DomLabels, Features] = GenDom_20Newsgroup(g_i, ClusSize, ComRate, NewsgroupIDs, ComIDs)

%%% Input
%
% g_i: the number of domain-specific networks within the same main cluster
% ClusSize: the domain cluster sizes
% ComRate: the common node ratio between domain-specific networks
% NewsgroupIDs: the newsgroup IDs of one category
% ComIDs: the common node IDs

%% Load dataset

load('../ExampleDataset/20Newsgroup.mat');

%% Extract documents from the newsgroups

NewsgroupPool = cell(length(NewsgroupIDs),1);

for i = 1:length(NewsgroupIDs)
    
    NewsgroupPool{i} = find(label == NewsgroupIDs(i));
    
end

%% Sample documents for domain-specific networks

DomIDs = cell(g_i,1);

for i = 1:g_i
    
    DomSample = [];
    
    for j = 1:length(NewsgroupPool)
        
        Idx = randsample(NewsgroupPool{j}, ClusSize);
        NewsgroupPool{j} = setdiff(NewsgroupPool{j}, Idx);
        DomSample = [DomSample; Idx];
        
    end
    
    DomIDs{i} = DomSample;
    
end

%% Construct domain-specific networks

DomNets = cell(g_i,1);
DomLabels = cell(g_i,1);
Features = cell(g_i,1);

for i = 1:g_i
    
    DomID_i = DomIDs{i};
    W_i = W(:,DomID_i');
    W_i_norms = zeros(size(W_i,2),1);
    
    for j = 1:size(W_i,2)
        
        W_i_norms(j) = norm(W_i(:,j),'fro');
        
    end
    
    Num = W_i_norms*W_i_norms';
    DomNet_i = ((W_i')*W_i)./Num;
    DomNet_i = DomNet_i - diag(diag(DomNet_i));
    DomLabel_i = label(DomID_i);
    
    % Remove empty rows and columns
    
    KeepRows= any(DomNet_i,2);
    DomNet_i = DomNet_i(KeepRows, KeepRows');
    DomID_i = DomID_i(KeepRows);
    DomLabel_i = DomLabel_i(KeepRows);
    W_i = W_i(:,KeepRows');
    
    DomNets{i} = DomNet_i;
    DomIDs{i} = DomID_i;
    DomLabels{i} = DomLabel_i;
    Features{i} = W_i;
    
end

% Generate common node IDs

for i = 1:g_i
    
    DomLabel_i = DomLabels{i};
    DomID_i = DomIDs{i};
    
    for j = 1:length(NewsgroupIDs)
        
        ComID_j = ComIDs{j};
        Idx = find(DomLabel_i == NewsgroupIDs(j));
        ComSize = round(ComRate*length(Idx));
        ComIdx = randsample(Idx, ComSize);
        DomID_i(ComIdx) = ComID_j(1:ComSize);
        
    end
    
    DomIDs{i} = DomID_i;
    
end

end