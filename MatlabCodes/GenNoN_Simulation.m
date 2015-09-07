%% Generate Synthetic NoN

function [DomNets, DomIDs, DomLabels, MainNet] = GenNoN_Simulation(mu, sigma)

%%% Input
%
% mu: the mean of random nodes adding/removing rate
% sigma: the standard deviation of random nodes adding/removing rate

%% Parameter Initialization

gs = [3,3,4]; % the main cluster sizes
ts = [5,5,5]; % the number of domain clusters

k = length(gs);
g = sum(gs);

DomNets = {};
DomIDs = {};
DomLabels = {};

IDs = cell(k,1); % The set of Common node IDs

for i = 1:k
    
    IDs{i} = 1:200;
    IDs{i} = IDs{i}';
    
end

NoiseIDs = cell(k,1); % The set of noisy domain node IDs

for i = 1:k
    
    NoiseIDs{i} = (i*1000 + 1):((i+1)*1000);
    NoiseIDs{i} = NoiseIDs{i}';
    
end

%% Generate a main network

MainNet = 0.5*ones(g,g);

for i = 1:k
    
    s = sum(gs(1:i-1))+1;
    e = sum(gs(1:i));
    MainNet(s:e,s:e) = 1;
    
end

MainNet = sparse(MainNet);

%% Generate the domain-specific networks

for i = 1:k
    
    [DomNet_i, DomID_i, DomLabel_i] = ...
        GenDom_Simulation(mu, sigma, gs(i), ts(i), IDs{i}, NoiseIDs{i});
    
    DomNets = [DomNets; DomNet_i];
    DomIDs = [DomIDs; DomID_i];
    DomLabels = [DomLabels; DomLabel_i];
    
end

end