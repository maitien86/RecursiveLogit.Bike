% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global maxstates;
    global isLinkSizeInclusive;
    % Get the LL value - How to get it:
    % ----------------------------------------------------
    % If Link size is included
    if (isLinkSizeInclusive == true)
        [LL, grad] = getODspecLL();
        return;
    end
    mu = 1; % MU IS NORMALIZED TO ONE
    % TO DO: compute this once and send these as parameters to this function
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    
    Mfull = getM(Op.x); % before it had two arguments (Op.x,false)
    [p q] = size(Mfull);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x); % before it had two arguments (Op.x,false)
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    % Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    % Initialize
    expV0 = zeros(1, q + 1);
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    gradV0 = zeros(q, Op.n);
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % LOOp ON ALL OBSERVATIONS
    % Compute the LL and gradient.
    for n = 1:nbobs       
        dest = Obs(n, 1);
        if true %(expV0(dest) == 0)   
            M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,dest);% addColumn;
            [expV, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters not fesible')
                return; 
            end                      
            expV0(dest) = expV(Obs(n,2)); % Get the exp(V0)
            gradV0(dest,:) = getGradV0(M, AttLc, Op, expV, Obs(n,2));
        else
        end
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        seq = 2;
        a = Obs(n,seq+1); % action state after origin
        % Value function at origin
        lnPn = - 1 * ((1/mu) * log(expV0(dest)));
        Gradient(n,:) = - gradV0(dest,:);
        while ((a ~= 0) && (seq < (maxstates - 1)))
            k = Obs(n,seq); % current state
            if (a > (lastIndexNetworkState+1))
                a = lastIndexNetworkState + 1; % The dest index for computation is always the same
            end        
            sumInstU = sumInstU + U(k,a);
            for i = 1:Op.n
                sumInstX(i) = sumInstX(i) + Atts(i).Value(k,a);
                %Un test
                fprintf('\nTest attribute value= %d \n', Atts(i).Value(k,a))
            end          
            seq = seq + 1;
            a = Obs(n,seq+1);
        end
        Gradient(n,:) = - gradV0(dest,:) + sumInstX;
        lnPn = lnPn + ((1/mu)*sumInstU) ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

%%
% Compute loglikelihood value with link size attribute 
%-----------------------------------------------------
function [LL, grad] = getODspecLL()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global maxstates;
    global LSatt;
    global LinkSize;
    
    % Get the LL value
    % ----------------------------------------------------
  
    mu = 1; % MU IS NORMALIZED TO ONE
    % TO DO: compute this once and send these as parameters to this function
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    [p q] = size(incidenceFull);
    LL = 0;
    grad = zeros(1, Op.n);
    % For the OD independence attributes
    for i = 1 : Op.n - 1
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % LoOp over all observation
    for n = 1:nbobs
        dest = Obs(n, 1);
        orig = Obs(n, 2);
        % get Link Size attributes
        LinkSize = LSatt(n).value;      
        AttLc(Op.n) =  Matrix2D(LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(Op.n).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(Op.n).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));       
        if true   
            % Get M and U
            Mfull = getM(Op.x, true); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            [Z, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters are not fesible')
                return; 
            end    
            expV0 = Z(orig);
            gradV0 = getGradV0(M, AttLc, Op, Z, orig);            
        end                
        % Get Utility
        Ufull = getU(Op.x, true); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);       
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        seq = 2;
        a = Obs(n,seq+1); % action state after origin
        lnPn = - 1 * ((1/mu) * log(expV0));
        Gradient(n,:) = - gradV0;
        while ((a ~= 0) && (seq < (maxstates - 1)))
            k = Obs(n,seq); % current state
            if (a > (lastIndexNetworkState+1))
                a = lastIndexNetworkState + 1; % The dest index for computation is always the same
            end        
            sumInstU = sumInstU + U(k,a);          
            for i = 1 : Op.n
                sumInstX(i) = sumInstX(i) + AttLc(i).Value(k,a);
            end                 
            seq = seq + 1;
            a = Obs(n,seq+1);
        end
        Gradient(n,:) = - gradV0 + sumInstX;
        lnPn = lnPn + ((1/mu)*sumInstU) ;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end