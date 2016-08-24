function [nodeResid, varParaOut, BICout, NWvarTopo] = gL0LS_SC_VAR...
    (Y, ar_o_max, h_initSearchScale, h_numIter, L1_init, ...
    ignoreSelfNodePen, EBIC_gamma, maxIter, cvg_tol, vb, DO_PARALLEL)


if DO_PARALLEL
    vb=0; %no verbose if parallel proc.
%    matlabpool(DO_PARALLEL)
    pool=parpool(DO_PARALLEL);
end

Ym = Y((ar_o_max+1):end,:);
[T,nNodes] = size(Y);
Tm = T - ar_o_max;
ar_rng = 1:ar_o_max;

hInit_L0 = 2*log(nNodes);
h_rng_L0 = logspace(log10(hInit_L0/h_initSearchScale), log10(hInit_L0*h_initSearchScale), h_numIter);
hInit_L1 = sqrt(2*log(nNodes));
h_rng_L1 = logspace(log10(hInit_L1/h_initSearchScale), log10(hInit_L1*h_initSearchScale), h_numIter);

H = setupNodeVarCell(Y, ar_o_max);

nodeResid = NaN*zeros(Tm,nNodes);
varParaOut = cell(nNodes,1);
BICout = cell(nNodes,1);
NWvarTopo = NaN*zeros(nNodes);
% for lp = 1:nNodes
parfor lp = 1:nNodes
    ignorePenaltyIdx = zeros(nNodes,1); 
    if ignoreSelfNodePen
        ignorePenaltyIdx(lp) = 1;
    end
    
    H_buf = H; % parfor warning
    
    BIC = struct;
    BIC.calc = zeros(ar_o_max, h_numIter);
    BIC_L1 = struct;
    BIC_L1.calc = [];
    for ar_o = ar_rng
        X = zeros(Tm, nNodes*ar_o);
        elt = cell(nNodes, 1);
        for v=1:nNodes
            [~, Col_j] = factorCol(X(1,:), v, ar_o*ones(numel(elt),1));
            elt{v} = Col_j;
            X(:,elt{v}) = H_buf{v}(:,1:ar_o);
        end
        
        % Do X colGrpOrth here, once!!
        [Xorth, Xwts] = colGrpOrth(X, elt); % do X mc in here.? no, too expensive?
        meanX = mean(Xorth);
        Xorth = bsxfun(@minus, Xorth, meanX);

        
        y_buf = Ym(:, lp);
        y = bsxfun(@minus, y_buf, mean(y_buf));
        y = bsxfun(@times, y, 1./(sqrt(sum(abs(y).^2,1))));
        
        sig2 = cov(y(:)); % for gL0LS tuning
        
        assert(all(size(Xorth) == [T-ar_o_max, (nNodes)*ar_o]));
        assert(all(size(y) == [T-ar_o_max, 1]));
        
        [bInit, ~, ~] = tikhLS_sure(0, y, Xorth, [-3,3], 1e-2, 15);
        
        %%
        if L1_init
            hLim = [0.7 0.7]; % min max, log10scaling for iterations
            hSearchIterMaxIter = 10;
            BIC_L1 = gPenLS_hSearchNested...
                ('LONE', ar_o, BIC_L1, y, Xorth, Xwts, h_rng_L1, hLim, sig2, bInit, elt, ...
                ignorePenaltyIdx, EBIC_gamma, ...
                maxIter, hSearchIterMaxIter, cvg_tol, vb);
            bInit = BIC_L1.localSoln.B;
        end
        %%
        hLim = [0.7 0.7]; % min max, log10scaling for iterations
        hSearchIterMaxIter = 10;
        BIC = gPenLS_hSearchNested...
            ('LZERO', ar_o, BIC, y, Xorth, Xwts, h_rng_L0, hLim, sig2, bInit, elt, ...
            ignorePenaltyIdx, EBIC_gamma, ...
            maxIter, hSearchIterMaxIter, cvg_tol, vb);
        %%
    end

    BICout{lp} = BIC;
    varParaOut{lp} = reshape(BIC.soln.B_unnorm, numel(BIC.solnElt{1}), []);
    NWvarTopo(lp,:) = sum(full(varParaOut{lp}~=0), 1);
    nodeResid(:,lp) = BIC.solnInfo.resid;
end

if DO_PARALLEL
%    if matlabpool('size') > 0
%        matlabpool close,
%    end
    if pool.NumWorkers > 0
        delete(pool),
    end
end


end

function H = setupNodeVarCell(Y, ar_o_max)

[T,nNodes ]= size(Y);

H = cell(1,nNodes);
for lp = 1:nNodes
    sig = Y(:,lp);
    hbuf = NaN*zeros(T-ar_o_max,ar_o_max);
    for ar_lp = 1:ar_o_max
        hbuf(:,ar_lp) = sig((ar_o_max-ar_lp)+1:(end-ar_lp));
    end
    H{lp} = bsxfun(@minus, hbuf, mean(hbuf));
    H{lp} = bsxfun(@times, H{lp}, 1./(sqrt(sum(abs(H{lp}).^2,1))));
end
end


function BIC = gPenLS_hSearchNested...
    (sparseMethod, ar_o, BIC, y, Xorth, Xwts, h_rng_init, hLim, sig2, bInit, elt,...
    ignorePenaltyIdx, EBIC_gamma, maxIter, hSearchIterMaxIter, cvg_tol, vb)
% NOTE: BIC.foobar(ar_o, hIdx) does NOT necessarily have same values of h
% across different ar_o.

h_rng= h_rng_init;

if numel(hLim)==1 % preset.
    h_rng = hLim;
elseif isempty(hLim) || numel(hLim) > 2, % we didn't specify range.
    %logscale, same order of magnitude as h_rng_init limits
    hLimScalingMin = myrange([log10(h_rng(1)), log10(h_rng(end))]) / 2 ; 
    %logscale, same order of magnitude as h_rng_init limits
    hLimScalingMax = myrange([log10(h_rng(1)), log10(h_rng(end))]) / 2 ;
else
    hLimScalingMin = hLim(1);
    hLimScalingMax = hLim(2);
end
% initial search
[BIC, hOptEst, stopFlag] = gPenLS_hSearch...
    (sparseMethod, ar_o, BIC, y, Xorth, Xwts, h_rng, sig2, bInit, elt, ...
    ignorePenaltyIdx, EBIC_gamma, maxIter, cvg_tol, vb);

hSearchItr = 1;
while hSearchItr < hSearchIterMaxIter
    if hOptEst ~= min(h_rng) && hOptEst ~= max(h_rng) && hSearchItr >= 3, break, end
    if stopFlag, break, end
    hLimMin = log10(hOptEst/(10^hLimScalingMin));
    hLimMax = log10(hOptEst*(10^hLimScalingMax));
    h_rng = logspace(hLimMin, hLimMax, length(h_rng));
    
    [BIC, hOptEst, stopFlag] = gPenLS_hSearch...
        (sparseMethod, ar_o, BIC, y, Xorth, Xwts, h_rng, sig2, bInit, elt, ...
        ignorePenaltyIdx, EBIC_gamma, maxIter, cvg_tol, vb);
    hSearchItr = hSearchItr +1;
end

end

function [BIC, hOptEst_constAR, stopFlag] = gPenLS_hSearch...
    (sparseMethod, ar_o, BIC, y, Xorth, Xwts, h_rng, sig2, bInit, elt, ...
    ignorePenaltyIdx, EBIC_gamma, maxIter, cvg_tol, vb)
Tm = size(y,1);


% BIC.hUsed(ar_o,:) = h_rng;
hIdx = 1;
for h = h_rng(:)'
    switch upper(sparseMethod)
        case 'LZERO'
            [soln, solnInfo] = gL0LS_SC...
                (y, Xorth, h, sig2, bInit, elt, ignorePenaltyIdx, maxIter, cvg_tol, Xwts, vb);
        case 'LONE'
            [soln, solnInfo] = gL1LS_SC...
                (y, Xorth, h, sig2, bInit, elt, ignorePenaltyIdx, maxIter, cvg_tol, Xwts, vb);
        otherwise
            error('Specify LZERO or LONE')
    end
    
    RSS = sum(solnInfo.resid.^2);
%     BIC.RSS(ar_o, hIdx) = sum(solnInfo.resid.^2);
%     BIC.T_effective(ar_o, hIdx) = Tm;
    gPenLSnorm = sum(solnInfo.gPenLSnorm_v);
%     BIC.gPenLSnorm(ar_o, hIdx) = sum(solnInfo.gPenLSnorm_v);
%     BIC.ar_ord(ar_o, hIdx) = ar_o;
        

    tauSj = size(Xorth,2); % or nchoosek(size(X,2), BIC.gPenLSnorm(ar_o, hIdx))
    BICbuf = ...
        Tm.* log(RSS./(Tm)) + ...
        gPenLSnorm.*ar_o.*( log(Tm) +  2*EBIC_gamma*log(tauSj));

    %overall best solution so far
    if (BICbuf < min(BIC.calc(BIC.calc~=0)))
        BIC.soln = soln;
        BIC.solnInfo = solnInfo;
        BIC.solnElt = elt;
    elseif isempty(BIC.calc)
        BIC.soln = soln;
        BIC.solnInfo = solnInfo;
        BIC.solnElt = elt;
    end

    % local solution, for giving (optional) initial L1 solution for 
    % subsequent L0 calc.
    BIC_opt_h = [];
    try BIC_opt_h = BIC.calc(ar_o,:); catch ,end 
    if BICbuf < min(BIC_opt_h(BIC_opt_h~=0))
        BIC.localSoln = soln;
    elseif isempty(BIC.calc)
        BIC.localSoln = soln;
    end
    
    BIC.calc(ar_o, hIdx) = BICbuf;
    
    hIdx = hIdx+1;
end
% for case when we have multiple points as minimum.
% find hOptEst_constAR preferentially as endpoints, otherwise doesn't matter.
bicbuf = BIC.calc(ar_o,:);
stopFlag = 0;
if all(min(bicbuf) == bicbuf), stopFlag = 1; hOptEst_constAR = h_rng(1); return;
elseif min(bicbuf) == bicbuf(1), 
    hOptEst_constAR = h_rng(1);
else
    hOptEst_constAR = h_rng( find(bicbuf == min(bicbuf), 1, 'last'));
end

end



function y = myrange(x)
    y = max(x) - min(x);
end

