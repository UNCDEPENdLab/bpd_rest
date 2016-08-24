
function [soln, solnInfo] = gL0LS_SC...
    (y, X, h_in, sig2, bInit, elt, ignorePenaltyIdx, maxIter, cvg_tol, Xwts, vb)

ignorePenaltyIdx = logical(ignorePenaltyIdx); % semantics, indicator variable

V = numel(elt);

if isempty(Xwts) % otherwise assume X is colGrpOrth already...
    [X, Xwts] = colGrpOrth(X, elt); % do X mc in here.? no, too expensive?
    meanX = mean(X);
    X = bsxfun(@minus, X, meanX);
end

meanY = mean(y);
y = bsxfun(@minus, y, meanY);

% optimisation
h = sig2*h_in;

[B, costFnOut, gL0norm_v] = gL0LS_calc...
    (y, X, bInit, h, V, elt, ignorePenaltyIdx, maxIter, cvg_tol, vb);

soln.B = sparse(B);
soln.B_unnorm = sparse(Xwts*B);

solnInfo.costFn = costFnOut;
solnInfo.gPenLSnorm_v = gL0norm_v;
% should we be judging e.g. BIC on scaled or unscaled model fit? because
% unscaled model is really what we're interested in.......?
solnInfo.resid = y - X*B;
% solnInfo.h_in = h_in;
end

%%%%% INTERNAL FUNCTIONS %%%%%%%
%%%%% INTERNAL FUNCTIONS %%%%%%%
%%%%% INTERNAL FUNCTIONS %%%%%%%
%%%%% INTERNAL FUNCTIONS %%%%%%%

function [B, costFnOut, gL0norm_v] = gL0LS_calc...
    (y, X, bInit, h, V, elt, ignorePenaltyIdx, maxIter, cvg_tol, vb)

B = bInit;
dotXB = X*B;
gL0norm_v = NaN*zeros(V,1);
for v=1:V
    gL0norm_v(v) = norm(B(elt{v}, :)) > 0;
end
activeSet = find(gL0norm_v);
iterSet = 1:V;
% iterSet = activeSet;
costFn = NaN*zeros(maxIter,1);

costFn(1) = costFnCalc(y, dotXB, gL0norm_v, h, ignorePenaltyIdx);

msg=[];
if vb,
    foo = struct(...
        'itr','Iter',...
        'costFn','CostFn',...
        'costFnSt','Cost Step',...
        'normVal','gLoNorm');
    msg = dispProgress(foo,msg);
    foo = struct(...
        'itr','init',...
        'costFn',num2str(costFn(1), 5),...
        'costFnSt','N/A',...
        'normVal',num2str(sum(gL0norm_v), 5));
    msg = dispProgress(foo,msg);
end

itr = 2;
while itr < maxIter
    for v = iterSet(:)'
        elt_v = elt{v};
        Xv = X(:,elt_v);
        Bv_prev = B(elt_v);
        err = y - dotXB;
        foo = Xv'*err + Bv_prev;
        if ignorePenaltyIdx(v)
            Bv_new = foo;
        else
            Bv_new = foo*(sqrt(foo'*foo) >= sqrt(h));
        end
%         dotXB = dotXB - Xv*Bv_prev + Xv*Bv_new;
        dotXB = dotXB -Xv*(Bv_prev-Bv_new);
        
        gL0norm_v(v) = norm(Bv_new) > 0; % boolean
        B(elt_v) = Bv_new;
    end
    
    costFn(itr) = costFnCalc(y, dotXB, gL0norm_v, h, ignorePenaltyIdx);
    
    
    if vb, 
        foo = struct(...
        'itr',num2str(itr),...
        'costFn',num2str(costFn(itr), 5),...
        'costFnSt',num2str(costFn(itr) - costFn(itr-1), 5),...
        'normVal',num2str(sum(gL0norm_v), 5));
        msg = dispProgress(foo,msg); 
    end
    
    [CVG, iterSet, activeSet] ...
        = cvgUpdate(itr, costFn, cvg_tol, gL0norm_v, iterSet, ...
        activeSet, V, B, X, maxIter, vb);
    if CVG, break, end
    if itr==maxIter,
        if vb, fprintf('\nIteration limit, did NOT converge\n\n'); end
        break
    end
    itr = itr+1;
end

costFnOut = costFn(itr);

end


function cFn = costFnCalc(y, dotXB, gL0norm_v, h, ignorePenaltyIdx)
cFn = (y-dotXB)'*(y-dotXB) + h*sum(gL0norm_v(ignorePenaltyIdx==0));
end

function [CVG, iterSet, activeSet, dotXB] = cvgUpdate(itr, costFn, cvg_tol, ...
    gL0norm_v, iterSet, activeSetPrev, V, B, X, maxIter, vb)
CVG = false;
activeSet = find(gL0norm_v);
if (abs(costFn(itr) - costFn(itr-1)) < cvg_tol)% if converged, locally or globally
    if (numel(iterSet) == V)
        if vb, fprintf('\n\n-- CONVERGED --\n\n'); end
        CVG = true;
    else % local converge, send back for global check
        iterSet = 1:V;
        dotXB = X*B; % update here to reduce numerical error
        if vb, fprintf('Local converge, check global\n'); end
    end
else % not converged
    if (numel(iterSet)== V) % just checked global
        if isequal(activeSet,activeSetPrev)
            % nothing changed in active set
            % so iterate over this as the iteration
            iterSet = find(gL0norm_v);
            if numel(iterSet) ~= V % in case we haven't been anywhere, be quiet!
                if vb, fprintf('Active set same after global search. Iter reduced set\n');end
            end
        else
            % passive set changed, so unstable. check global again.
            iterSet = 1:V;
            if vb, fprintf('Active set changed, check global\n');end
        end
    elseif itr +3 >=maxIter
        iterSet = 1:V;
        %panic, almost out of iterations... check global.
        dotXB = X*B; % update here to reduce numerical error
    else
        % iterate reduced set, no change
    end
end

end

%%%%% OTHER FUNCTIONS %%%%%%
%%%%% OTHER FUNCTIONS %%%%%%
%%%%% OTHER FUNCTIONS %%%%%%
%%%%% OTHER FUNCTIONS %%%%%%

function [An, AnWts] = colGrpOrth(A,elt)
Ngrp = numel(elt);
An = A; % should be same size if sparse
AnWts_j = cell(Ngrp,1);
for j = 1:Ngrp
    eltj = elt{j};
    Aj = A(:,eltj);
    cAjAj = chol(Aj.'*Aj);
    AnWts_j{j} = inv(cAjAj);
    An(:,eltj) = Aj/cAjAj;
end
AnWts = sparse(blkdiag(AnWts_j{:}));
end

function msg = dispProgress(in, msg)
msg = logProc(msg,[...
    sprintf('%4s',in.itr) ...
    sprintf('%15s', in.costFn) ...
    sprintf('%15s', in.costFnSt) ...
    sprintf('%9s', in.normVal) ...
    ]);
end

function [msg] = logProc(msg,proc)

% display in command window
disp(proc);

% store in msg cell array
% msg = {msg{:},proc};
msg = {msg {proc}};

end