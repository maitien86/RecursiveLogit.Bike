
function P = getP(expV, M)

% %     clear P;
%      for k=1:length(expV)
%         P(k,:) = (M(k,:) .* expV') / expV(k) ;
%         testRowsum = P(k,:) * ones(length(expV),1);
%      end
%      P = sparse(P);

    %---------- New idea Tien -------------
    N = size(M,1);
    MI = M; 
    MI(find(M)) = 1;
    
%     %some elements of invZdiag may be inf because too large:
% when computing link flows uncomment    
%    expV(expV<realmin)=realmin;
    invexpV=1./expV;
    
    Zdiag =  spdiags(expV,0,N,N);
    invZdiag =  spdiags(invexpV,0,N,N);%invexpV instead of 1./expV
    
    P = M .* (invZdiag * MI * Zdiag);
end

