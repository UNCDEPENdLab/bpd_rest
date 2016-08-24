function [nullThresh] = SCoTMI_H0calc(pval, nTimePoints, nNodes, paraVar, noiseVar, ar_sim_ord, nRep)


%defaults
if ~exist('paraVar', 'var') || isempty(paraVar),  paraVar = 0.05; end
if ~exist('noiseVar', 'var') || isempty(noiseVar),  noiseVar = 0.05; end
if ~exist('ar_sim_ord', 'var') || isempty(ar_sim_ord),  ar_sim_ord = 2; end
if ~exist('nRep', 'var') || isempty(nRep),  nRep = 50; end


opt.vb=0;
opt.ignoreSelfNodePen = 1;

opt.NFFT = 2^nextpow2(500);

opt.ar_o_max = 4;
% opt.EBIC_gamma = 0.1;

%%
nCmpr = nNodes*(nNodes-1)/2;

ctmi_cat = NaN*zeros(nCmpr, nRep);

NULL_UNIVARIATE_ARFIT = true;
for rep = 1:nRep
[ctMI_o, ~,~, ~, ~] = ...
     calc_H0...
     (nTimePoints, nNodes, paraVar, noiseVar, ar_sim_ord, opt, NULL_UNIVARIATE_ARFIT);

foo = ctMI_o(:);
foo(isnan(foo)) = [];
ctmi_cat(:,rep) = foo;

end    
%%

meanEst_perRep = mean(ctmi_cat);
stdErr_meanEst = std(meanEst_perRep);
meanEst = mean(ctmi_cat(:));

stdEst_perRep = std(ctmi_cat);
stdErr_stdEst = std(stdEst_perRep);
stdEst= std(ctmi_cat(:));

nullThresh = norminv(1-(pval/nCmpr), meanEst, stdEst);

end


function [ctMI, NWvarTopo, varPara, CHOLFLAG, nodeResid] = ...
    calc_H0...
        (nSamples, nNodes, paraVar, noiseVar, ar_sim_o, opt, NULL_UNIVARIATE_ARFIT)


Y = NaN*zeros(nSamples, nNodes);

B= 1;
A = [ones(nNodes,1) randn(nNodes,ar_sim_o)./sqrt(1./paraVar) ];
for lp = 1:nNodes
    lpCount = 1;
    while lpCount < 100
    
        if all(abs(roots(A(lp,:))) < 1)
            break
        else % try again...
%             disp(lp)
%             disp(abs(roots(A(lp,:))))
            A(lp,:) = [1 randn(1,ar_sim_o)./sqrt(1./paraVar)];
        end
        lpCount = lpCount+1;
    end
    
    initSkip = ar_sim_o; %number to skip, to remove startup artefacts 
    Ybuf = filter(B, A(lp,:), randn(initSkip + nSamples,1)./sqrt(1./noiseVar));
    Y(:,lp) = Ybuf(end-nSamples+1:end);
end

if NULL_UNIVARIATE_ARFIT
    [ctMI, NWvarTopo, varPara, BIC, CHOLFLAG, ...
        nodeResid, partCoh, Aomega] ...
        = SCoTMI_threshold_tmp(Y, opt);
else
    warning('Use univariate model for the null hypthesis test')
%     [ctMI, NWvarTopo, varPara, BIC, CHOLFLAG, nodeResid] = ctMI_calc(Y, opt);
end

end