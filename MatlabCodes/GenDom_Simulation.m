%% Generate synthetic domain-specific networks within the same main cluster

function [DomNets, DomIDs, DomLabels] = GenDom_Simulation(mu, sigma, g_i, t, IDs, NoiseIDs)

%%% Input
%
% mu: the mean of random nodes adding/removing rate
% sigma: the standard deviation of random nodes adding/removing rate
% g_i: the number of domain-specific networks within the same main cluster
% t: the number of domain clusters
% IDs: the IDs of common nodes
% NoiseIDs: the IDs of noisy domain nodes

%% Generate the underlying clustering structure

N = length(IDs);
TrueID = zeros(N,1);
TrueLabel = zeros(N,1);
u = []; % the end nodes of edges
v = []; % the end nodes of edges

for i = 1:t
    
    Clus_i_size = round(N/t);
    s = (i-1)*Clus_i_size+1;
    e = i*Clus_i_size;
    Clus_i_idx = s:e;
    
    u = [u; kron(Clus_i_idx', ones(Clus_i_size,1))];
    v = [v; kron(ones(Clus_i_size,1), Clus_i_idx')];
    
    SampleIDs = randsample(IDs,Clus_i_size);
    TrueID(s:e) = SampleIDs;
    IDs = setdiff(IDs,SampleIDs);
    TrueLabel(s:e) = i;
    
end

TrueNet = sparse(u, v, ones(length(u),1), N, N);
TrueNet = TrueNet - diag(diag(TrueNet));

%% Generate synthetic domain-specific networks

FlipRate0 = 0.8;
FlipRate1 = 0.05;

DomNets = cell(g_i,1);
DomIDs = cell(g_i,1);
DomLabels = cell(g_i,1);

for i = 1:g_i
    
    % Randomly flip 1 to 0 and flip 0 to 1
    
    rk = 0;
    
    while rk ~= N
        
        InputNet = TrueNet;
        InputID = TrueID;
        InputLabel = TrueLabel;
        
        ExistEdge = find(InputNet);
        NonExistEdge = find(~InputNet);
        
        % Randomly flip 1 to 0
        
        FlipEdge = randsample(ExistEdge,round(FlipRate0*length(ExistEdge)));
        InputNet(FlipEdge) = 0;
        
        % Randomly flip 0 to 1
        
        FlipEdge = randsample(NonExistEdge,round(FlipRate1*length(NonExistEdge)));
        InputNet(FlipEdge) = 1;
        
        % Make the matrix symmetric
        
        InputNet = triu(InputNet) + triu(InputNet)' - 2*diag(diag(InputNet));
        InputNet(InputNet>0) = 1;
        
        rk = length(find(any(InputNet,2)));
        
    end
    
    % Randomly remove and add nodes
    
    DomNets{i} = InputNet;
    DomIDs{i} = InputID;
    DomLabels{i} = InputLabel;
    
    % Randomly remove nodes

    RemoveRate = 0; % RemoveRate follows normal distribution with mean mu and standard deviation sigma
    
    if mu ~= 0 || sigma ~= 0
        
        while RemoveRate <= 0 || RemoveRate >= 1 || round(RemoveRate*N) == N
            
            RemoveRate = mu + sigma*randn(1);
            
        end
        
    end
    
    RemoveNodes = randsample(1:N,round(RemoveRate*N));
    DomNets{i}(RemoveNodes',:) = [];
    DomNets{i}(:,RemoveNodes) = [];
    DomIDs{i}(RemoveNodes) = [];
    DomLabels{i}(RemoveNodes) = [];
    
    % Randomly add nodes
    
    AddRate = 0; % AddRate follows normal distribution with mean mu and standard deviation sigma
    
    if mu ~= 0 || sigma ~= 0
        
        while AddRate <= 0 || AddRate >= 1 || round(AddRate*N) == N
            
            AddRate = mu + sigma*randn(1);
            
        end
        
    end
    
    AddNodes = round(N*AddRate);
    DomNets{i}(end+1:end+AddNodes,:) = 0;
    DomNets{i}(:,end+1:end+AddNodes) = 0;
    tmpID = randsample(NoiseIDs,AddNodes);
    DomIDs{i}(end+1:end+AddNodes) = tmpID;
    NoiseIDs = setdiff(NoiseIDs,tmpID);
    DomLabels{i}(end+1:end+AddNodes) = 0;
    
    rk = 0;
    
    while rk ~= size(DomNets{i},1)
        
        TmpBlock = sparse(size(DomNets{i},1),AddNodes);
        FlipEdge = randsample(find(~TmpBlock),round(FlipRate1*length(find(~TmpBlock))));
        TmpBlock(FlipEdge) = 1;
        DomNets{i}(:,end-AddNodes+1:end) = TmpBlock;
        
        DomNets{i} = triu(DomNets{i}) + triu(DomNets{i})' - 2*diag(diag(DomNets{i}));
        DomNets{i}(DomNets{i}>0) = 1;
        
        rk = length(find(any(DomNets{i},2)));
        
    end
    
end

end