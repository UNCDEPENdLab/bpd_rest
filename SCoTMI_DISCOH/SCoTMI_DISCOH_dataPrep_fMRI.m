% Extracts the first principal component of the voxel time series 
% from each node. credit to svdsecon (Vipin Vijayan (2014))
%
% Inputs:
%       fullData : time x voxels
%       nodeMap  : voxels x 1, mapping which voxels to group into nodes
%       u_nodes  : the actual nodes to use (may be a subset from nodeMap)
% 
% Output:
%       Y        : time x nodes
%
function [Y] = SCoTMI_DISCOH_dataPrep_fMRI(fullData, nodeMap, u_nodes)
T = size(fullData,1);
NNodes = length(u_nodes);

Y = zeros(T,NNodes);
for lp = 1:NNodes
    nodeData = fullData(:,nodeMap(:)==u_nodes(lp));
    X = bsxfun(@minus, nodeData, mean(nodeData));
    [m,n] = size(X);
    
    if  m <= n
        C = X*X';
        [U,~] = eigs(C,1);
    else
        C = X'*X;
        [V,D] = eigs(C,1);
        U = X*V;
        s = sqrt(abs(diag(D)));
        U = bsxfun(@(x,c)x./c, U, s');
    end
    Y(:,lp) = U; 
end