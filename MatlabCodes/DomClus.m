%% Simultaneous clustering of the domain-specific networks

function Us = DomClus(As, H, D, O, n_v, a, t_u, t_v, MaxIter, epsilon)

%%% Input
%
% As: a set of domain-specific networks
% H: the main cluster indicator matrix
% D: a mapping matrix
% O: a mapping matrix
% n_v: a vector of size k of the number of nodes in V
% a: a regularization parameter of NoNClus
% t_u: a vector of the numbers of domain clusters in domain-specific networks
% t_v: a vector of the numbers of domain clusters in hidden factor matrices
% MaxIter: the maximal number of iterations for alternating minimization
% epsilon: the convergence parameter

%% Normalize the domain specific networks

g = size(H,1);
k = size(H,2);
n_u = zeros(g,1);

for i = 1:g
    
    As{i} = As{i}/sqrt(trace(As{i}'*As{i}));
    n_u(i) = size(As{i},1);
    
end

%% Initialization

Us = cell(g,1);
Vs = cell(k,1);

for i = 1:g
    
    Us{i} = rand(n_u(i),t_u(i));
    Us{i} = Us{i}/sqrt(trace(Us{i}'*Us{i}));
    
end

for i = 1:k
    
    Vs{i} = rand(n_v(i),t_v(i));
    Vs{i} = Vs{i}/sqrt(trace(Vs{i}'*Vs{i}));
    
end

J1 = obj(As, Us, Vs, H, D, O, g, k, a);
Iter = 1;
Delta = 99999;
Objs = J1;
Deltas = Delta;

%% Alternating update

while Delta > epsilon && Iter <= MaxIter
    
    % Update U
    
    for i = 1:g
        
        Num = As{i}*Us{i};
        Denom = (Us{i})'*Us{i};
        Denom = Us{i}*Denom;
        
        for j = 1:k
            
            % Compute W
            
            W = sparse(D{i,j}*Us{i});
            W = (O{i,j}')*W;
            W = sparse((Vs{j}')*W);
            W = sparse(Vs{j}*W);
            W = O{i,j}*W;
            W = (D{i,j}')*W;
            Num = Num + a*H(i,j)*W;
            
            % Compute Y
            
            Y = sparse(D{i,j}*Us{i});
            Y = (D{i,j}')*Y;
            Y = sparse((Us{i}')*Y);
            Y = sparse(Us{i}*Y);
            Y = D{i,j}*Y;
            Y = (D{i,j}')*Y;
            Denom = Denom + a*H(i,j)*Y;
            
        end
        
        Us{i} = Us{i}.*((Num./Denom).^(0.25));
        
    end
    
    % Update V
    
    for i = 1:k
        
        Num = sparse(n_v(i),t_v(i));
        Denom = sparse(n_v(i),t_v(i));
        
        for j = 1:g
            
            % Compute Q
            
            TmpMat = sparse(D{j,i}*Us{j});
            TmpMat = (O{j,i}')*TmpMat;
            Q = sparse((TmpMat')*Vs{i});
            Q = TmpMat*Q;
            Num = Num + H(j,i)*Q;
            
            % Compute R
            
            TmpMat = sparse(O{j,i}*Vs{i});
            TmpMat = (O{j,i}')*TmpMat;
            R = sparse((TmpMat')*Vs{i});
            R = TmpMat*R;
            Denom = Denom + H(j,i)*R;
            
        end
        
        Vs{i} = Vs{i}.*((Num./Denom).^(0.25));
        
    end
    
    J2 = obj(As, Us, Vs, H, D, O, g, k, a);
    Delta = J1 - J2;
    Objs = [Objs, J2];
    Deltas = [Deltas, Delta];
    J1 = J2;
    Iter = Iter + 1;
    
end

end

%% Objective function

function J = obj(As, Us, Vs, H, D, O, g, k, a)

J = 0;

for i = 1:g
    
    J = J + norm(As{i} - Us{i}*Us{i}', 'fro')^2;
    
    for j = i:k
        
        ProjU = D{i,j}*Us{i};
        ProjV = O{i,j}*Vs{j};
        J = J + a*H(i,j)*(norm(ProjU*(ProjU)' - ProjV*(ProjV)', 'fro')^2);
        
    end
    
end

end