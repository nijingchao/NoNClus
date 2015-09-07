%% Generate 20newsgroup NoN

function [DomNets, DomIDs, DomLabels, MainNet] = GenNoN_20Newsgroup(ComRate)

%%% Input
%
% ComRate: the common node ratio between domain-specific networks

%% Parameter Initialization

gs = [10,10,10]; % the main cluster sizes
ClusSize = 50; % the domain cluster sizes
Newsgroups = {[2,3,4,5], [8,9,10,11], [17,18,19,20]}; % The newsgroup IDs of three categories

k = length(gs);
g = sum(gs);

DomNets = {};
DomIDs = {};
DomLabels = {};
Features = {};

ComIDs = {[20001:1:20000+round(ComRate*ClusSize)]'; ...
    [22001:1:22000+round(ComRate*ClusSize)]'; ...
    [23001:1:23000+round(ComRate*ClusSize)]'; ...
    [24001:1:24000+round(ComRate*ClusSize)]'}; % The common node IDs

%% Generate domain-specific networks

for i = 1:k
    
    [DomNet_i, DomID_i, DomLabel_i, Feature_i] = ...
        GenDom_20Newsgroup(gs(i), ClusSize, ComRate, Newsgroups{i}, ComIDs);
    
    DomNets = [DomNets; DomNet_i];
    DomIDs = [DomIDs; DomID_i];
    DomLabels = [DomLabels; DomLabel_i];
    Features = [Features; Feature_i];
    
    % Shuffle common IDs
    
    AllComIDs = cell2mat(ComIDs);
    idx = randperm(length(AllComIDs));
    AllComIDs = AllComIDs(idx);
    ComIDs = mat2cell(AllComIDs, round(ComRate*ClusSize)*ones(1,4), 1);
    
end

%% Generate the main network

Main_norms = zeros(g,1);
MainFeatures = zeros(size(Features{1},1),g);

for i = 1:g
    
    MainFeature_i = sum(Features{i},2);
    Main_norms(i) = norm(MainFeature_i,'fro');
    MainFeatures(:,i) = MainFeature_i;
    
end

Num = Main_norms*Main_norms';
MainNet = ((MainFeatures')*MainFeatures)./Num;
MainNet = MainNet - diag(diag(MainNet));
MainNet = sparse(MainNet);

end