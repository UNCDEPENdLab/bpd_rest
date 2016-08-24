function [An, AnWts] = colGrpOrth(A,elt)
Ngrp = numel(elt);
An = zeros(size(A));
AnWts_j = cell(Ngrp,1);
for j = 1:Ngrp
    Aj = A(:,elt{j});
    cAjAj = chol(Aj.'*Aj);
    AnWts_j{j} = inv(cAjAj);
    An(:,elt{j}) = Aj/cAjAj;
end
AnWts = sparse(blkdiag(AnWts_j{:}));
end