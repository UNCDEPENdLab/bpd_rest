% if hLim is scalar, then we only calculate that. 
% if hLim is a 2-vector, we calculate logscale within that range, maxIter
% times.
% else use default values.

% TODO: 
% 
% change hLim, somehow to require normalised model and
% therefore more sensible values
%
% condLimit as optional process
%
% ?default for just tikh soln with specified reg-para., no SURE?

function [b, h_min, SURE] = tikhLS_sure(vb, y, G, hLim, condLimit, maxIter)
SURE = struct;
if isempty(maxIter), maxIter = 100; end
if isempty(condLimit), condLimit = 1e-2; end
if numel(hLim)==1 % preset.
    h_min = hLim;
elseif isempty(hLim) || numel(hLim) > 2, % we didn't specify range.
    hLimMin = -10; %logscale
    hLimMax = 10; %logscale
else
    hLimMin = hLim(1);
    hLimMax = hLim(2);
end

[U, S, V] = svd(G, 'econ');

[N,~] = size(G);
nDim = sum(diag(S)./S(1) > condLimit);
if vb, figure(); plot(diag(S)); title({['TIKH: Singular Values of weighted Gain'],['Cond:' num2str(cond(S),3) ', nDim:' num2str(nDim) '/' num2str(N)]}); end

SS = diag(S).^2;

if ~exist('h_min', 'var') || isempty(h_min); % check a range of values... otherwise use input
    vb_buf = vb;
    [h_min_init, SURE.init] = SUREest(hLimMin,hLimMax,maxIter, U, S, V, G, N, nDim, SS, vb_buf, y);
    hLimMin = log10(h_min_init/10);
    hLimMax = log10(h_min_init*10);
    vb_buf = vb;
    [h_min, SURE.exact] = SUREest(hLimMin,hLimMax,maxIter, U, S, V, G, N, nDim, SS, vb_buf, y);
end
% final solution
D = diag(S)./(SS + h_min);
if nDim < numel(diag(S)),
    D(nDim+1:end) = 0;
end
b = V*diag(D)*U.'*y;

    
end


function [h_min, SURE] = SUREest(hLimMin,hLimMax,maxIter, U, S, V, G, N, nDim, SS, vb_buf, y)
        SURE = NaN*zeros(maxIter, 2);
        iter = 1;
        for ii = logspace(hLimMin, hLimMax, maxIter)
            h = ii;
            D = diag(S)./(SS + h);
            if nDim < numel(diag(S)),
                D(nDim+1:end) = 0;
            end
            b = V*diag(D)*U.'*y;
            norm2res = norm(y - G*b, 2)^2;
            sigma2 = norm2res./N;
            SURE_ii = norm2res + 2*sigma2*sum(SS./(SS + h));
            SURE(iter,:) = [SURE_ii, h];
            iter = iter+1;
        end
        
        [~, minInd] = min(SURE(:,1));
        h_min = SURE(minInd, 2);
        if vb_buf,figure(); 
            loglog(SURE(:,2), SURE(:,1)), title(['Min h-tikh = ' num2str(h_min)]); drawnow; ...
            
            fprintf('Min h-tikh = %d\n',h_min) 
        end
    end