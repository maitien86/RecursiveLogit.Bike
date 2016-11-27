%   Compute Link Size attribute from data
%   
%%
function ExpV_is_ok = getLinkSizeAtt_new()
    global incidenceFull; 
    global Mfull;
    global Obs;     % Observation
    global nbobs;  
    global LSatt;
    global Atts;
    global EstimatedTime;
    global Op;
    %beta=[-1.94,-0.93,-0.85,-3.8,-1.6,0.33,1.49,0.7,-5.1,4.07,-0.22,0.34,-0.31,1.37];
    %beta = [-1,-1,-1, 0]';
    beta = [-5.2570;-4.6702;-5.0949;-1.5374;-5.1486]';
    %beta = [-1.0,0,0,0]';
    % ----------------------------------------------------
    mu = 1; % MU IS NORMALIZED TO ONE
    
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    %Mfull = getM(beta,false);
    
    u = beta(1) * Atts(1).value;
    for i = 2:5%Op.n
        u = u + beta(i) * Atts(i).value; 
    end
    expM = u;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    Mfull = expM;
    
    
    [p q] = size(Mfull);
    M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState); 
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
   
    for n = 1:nbobs   
        n
        dest = Obs(n, 1);
        orig = Obs(n, 2);
        if true %expV0(dest) == 0   
            M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,dest);
            [expV, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                ExpV_is_ok = false;
                disp('ExpV is not feasible')
                return; 
            end  
            P = getP(expV, M);   
        end
        G = sparse(zeros(size(expV)));
        G(orig) = 1;
        I = speye(size(P));
        F = (I-P')\G;                        
        if (min(F) < 0)                
            ToZero = find(F <= 0);
            for i=1:size(ToZero,1)
                F(ToZero(i)) = 1e-9;
            end
        end
        ips = F;
        ips = ips(1:lastIndexNetworkState); % dummy link should not be included
        
        %% Add correction term
%         TT = max(EstimatedTime); % vector
%         ips = ips(1:lastIndexNetworkState); 
%         TT =  TT(1:lastIndexNetworkState);       
%         ips = ips .* TT';
%         
        
        ips(size(incidenceFull,2),1) = 0;
        %lsatt = sparse(zeros(size(Atts(1).Value)));      
        I = find(incidenceFull);
        [nbnonzero, c] = size(I);
        ind1 = zeros(nbnonzero,1);
        ind2 = zeros(nbnonzero,1);
        s = zeros(nbnonzero,1);
        for i = 1:nbnonzero
            [k a] = ind2sub(size(incidenceFull), I(i));
            ind1(i) = k;
            ind2(i) = a;
            s(i) =  ips(a);
        end    
        lsatt = sparse(ind1, ind2, s, size(incidenceFull,1), size(incidenceFull,2));        
        LSatt(n).value = lsatt;
    end
    ExpV_is_ok = true;
end
