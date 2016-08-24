function [partCoh, CHOLFLAG] = partCohCalc(nodeResid, Aomega, NFFT)
%%
nNodes = size(nodeResid,2);

covResid = cov(nodeResid);

if all(eig(covResid) > 0)
    CHOLFLAG = 1;
    L = chol(covResid, 'lower');
    
    partCoh_cmp= cell(NFFT,1);
    for wInd = 1:NFFT
        Aw = NaN*zeros(nNodes);
        for nodeLP = 1:nNodes
            Aw(:,nodeLP) = Aomega{nodeLP}(wInd,:);
        end
        
        partCoh_Hmatrix = L\Aw;
        buf = partCoh_Hmatrix'*partCoh_Hmatrix;
        partCoh_cmp{wInd} = buf;
    end
    
%     partCoh_Hmatrix = cell(NFFT,1);
%     for wInd = 1:NFFT
%         partCoh_Hmatrix{wInd} = NaN*zeros(nNodes);
%         Aw = NaN*zeros(nNodes);
%         for nodeLP = 1:nNodes
%             Aw(:,nodeLP) = Aomega{nodeLP}(wInd,:);
%         end
%         partCoh_Hmatrix{wInd} = L\Aw;
%     end
%     
%     partCoh_cmp= cell(NFFT,1);
%     for wInd = 1:NFFT
% %         if mod(wInd,50)==0, disp(wInd),end
%         partCoh_cmp{wInd} = NaN*zeros(nNodes);
%         for nodeLP = 1:nNodes
%             for refLP = 1:nNodes
%                 if nodeLP < refLP
%                     continue
%                 end
%                 
%                 partCoh_cmp{wInd}(nodeLP,refLP) =...
%                     partCoh_Hmatrix{wInd}(:,nodeLP)' * partCoh_Hmatrix{wInd}(:,refLP);
%                 
%             end
%         end
%     end
    
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
    %diagonal matrix to pre- and post-multiply
    normMatSp = sparse(1:nNodes,1:nNodes, 1./sqrt(diag(partCoh_cmp{wInd})));
    partCoh_BUF = -normMatSp*partCoh_cmp{wInd}*normMatSp;
    
    partCoh_BUF(triu(partCoh_BUF,0)~=0) = NaN;
    partCoh(:,:,wInd) = partCoh_BUF;
end

% partCoh = cell(NFFT,1);
% for wInd = 1:NFFT
%     partCoh{wInd} = NaN*zeros(nNodes);
%     for nodeLP = 1:nNodes
%         for refLP = 1:nNodes
%             if nodeLP < refLP
%                 continue
%             end
%             
%             partCoh{wInd}(nodeLP,refLP) = ...
%                 -partCoh_cmp{wInd}(nodeLP,refLP) ...
%                     ./(sqrt(partCoh_cmp{wInd}(nodeLP,nodeLP)) *...
%                        sqrt(partCoh_cmp{wInd}(refLP,refLP)) ...
%                       );
%         end
%     end
% end

end