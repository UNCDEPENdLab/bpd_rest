function [ctMI, NWvarTopo, varPara, BIC, CHOLFLAG, nodeResid, partCoh, Aomega] ...
    = SCoTMI_threshold_tmp(Y, opt)


%DEFAULTS
if isempty(opt)
    fprintf('Using default optimisation parameters')
    
    opt.vb=1;

    opt.ignoreSelfNodePen = 1;
    
    opt.NFFT = 2^nextpow2(1e3);
    
    opt.ar_o_max = 5;

end

[T,nNodes] = size(Y);
Tm = T-opt.ar_o_max;

[nodeResid, varPara, BIC, NWvarTopo] ...
    = arNullNetworkModelFit(Y, opt.ar_o_max, opt.ignoreSelfNodePen, opt.vb);


%%
Aomega = cell(nNodes,1);
for nodeLP=1:nNodes
    % important:: fft explicitly along columns.
    Aomega{nodeLP} =    -   (1./sqrt(opt.NFFT))*fft(full(varPara{nodeLP}), opt.NFFT, 1);
    Aomega{nodeLP}(:,nodeLP) = Aomega{nodeLP}(:,nodeLP) + 1; %already did minus, now add 1.
end
%%
[partCoh, CHOLFLAG] = partCohCalc(nodeResid, Aomega, opt.NFFT);

ctMI = tMIcalc(partCoh, opt.NFFT, Tm);



end


function [partCoh, CHOLFLAG] = partCohCalc(nodeResid, Aomega, NFFT)
%%
nNodes = size(nodeResid,2);

covResid = cov(nodeResid);

if all(eig(covResid) > 0)
    CHOLFLAG = 1;
    L = chol(covResid, 'lower');
    
    partCoh_cmp= cell(NFFT,1);
%     partCoh_Hmatrix = cell(NFFT,1);
    for wInd = 1:NFFT
%         partCoh_Hmatrix{wInd} = NaN*zeros(nNodes);
        Aw = NaN*zeros(nNodes);
        for nodeLP = 1:nNodes
            Aw(:,nodeLP) = Aomega{nodeLP}(wInd,:);
        end
%         partCoh_Hmatrix{wInd} = L\Aw;
        
        partCoh_Hmatrix = L\Aw;
        buf = partCoh_Hmatrix'*partCoh_Hmatrix;
%         buf(triu(buf,1)~=0) = NaN;
        partCoh_cmp{wInd} = buf;
    end
    
else % Do the same explicitly (without Cholesky stabilisation and speedup).
    CHOLFLAG = 0;
    invCovResid = eye(nNodes)\covResid;
    
    partCoh_g = cell(NFFT,1);
    for wInd = 1:NFFT
        partCoh_g{wInd} = NaN*zeros(nNodes);
        Aw = NaN*zeros(nNodes);
        for nodeLP = 1:nNodes
            Aw(:,nodeLP) = Aomega{nodeLP}(wInd,:);
        end
        partCoh_g{wInd} = Aw;
    end
    
    
    partCoh_cmp= cell(NFFT,1);
    for wInd = 1:NFFT
        if mod(wInd,50)==0, disp(wInd),end
        partCoh_cmp{wInd} = NaN*zeros(nNodes);
        for nodeLP = 1:nNodes
            for refLP = 1:nNodes
                if nodeLP < refLP
                    continue
                end
                
                partCoh_cmp{wInd}(nodeLP,refLP) =...
                    partCoh_g{wInd}(:,nodeLP)' * invCovResid * partCoh_g{wInd}(:,refLP);
                
            end
        end
    end
end

partCoh = NaN*zeros(nNodes, nNodes, NFFT);
for wInd = 1:NFFT
%     normMat = diag(1./sqrt(diag(partCoh_cmp{wInd})));
%     partCoh_BUF = -normMat*partCoh_cmp{wInd}*normMat;
    normMatSp = sparse(1:nNodes,1:nNodes, 1./sqrt(diag(partCoh_cmp{wInd})));
    partCoh_BUF = -normMatSp*partCoh_cmp{wInd}*normMatSp;
    
    partCoh_BUF(triu(partCoh_BUF,0)~=0) = NaN;
    partCoh(:,:,wInd) = partCoh_BUF;
end



end

function tMI = tMIcalc(partCoh, NFFT, Tm)
%%

% nNodes = size(partCoh,1);
tMIscaling = -Tm./(2*NFFT);
tMIbuf = tMIscaling * log(1 - abs(partCoh).^2);
tMI = sum(tMIbuf, 3);


end

function [nodeResid, varPara, BICout, NWvarTopo] ...
    = arNullNetworkModelFit(Y, ar_o_max, ignoreSelfNodePen, vb)

[T,nNodes] = size(Y);
Tm = T - ar_o_max;

nodeResid = NaN*zeros(Tm,nNodes);
varPara = cell(nNodes,1);
BICout = cell(nNodes,1);
NWvarTopo = zeros(nNodes, nNodes);

minOrd = double(ignoreSelfNodePen );
maxOrd = ar_o_max;
for lp = 1:nNodes
    [resOut, phi_out, order_out, BIC] = ar_est_whitenoise(Y(:,lp), minOrd, maxOrd);
    nodeResid(:,lp) = resOut;
%     nodeResid(:,lp) = resOut(ar_o_max+1:end);
    if order_out >0
        varPara{lp} = sparse(1:order_out, lp, phi_out(2:end), order_out, nNodes);
    else
        varPara{lp} = sparse(1,1,0,1,nNodes);
    end
    BICout{lp} = BIC;
    NWvarTopo(lp,lp) = order_out;
end

% output varPara as a cell per node, ar_o*nNodes. but sparse, only 1col
% each.
end

function [resOut, phi_out, order_out, BIC] = ar_est_whitenoise(y, minOrd, maxOrd)
y = bsxfun(@minus, y(:), mean(y(:)));
phi_out = [1 0];
order = 0;
order_out= 0;
BIC = NaN*zeros(maxOrd+1, 1);
n = length(y);
if (order >= minOrd)
    BIC(1) = n*log(var(y));
    BICmin = BIC(1);
end


while order < maxOrd
    order = order+1;
    [phi, ~] = arburg(y, order);
    res = filter(phi(:), 1, y);
    res_trunc = res(maxOrd+1:end);
    nVar = (res_trunc'*res_trunc)/(n-maxOrd);
    BIC(order+1) = n*log(nVar) + order*log(n);
    
    if ~exist('BICmin', 'var') || ((BIC(order+1) < BICmin) && (order >= minOrd))
        resOut = res_trunc;
        BICmin = BIC(order+1);
        order_out = order;
        phi_out = phi;
    end

end

end

