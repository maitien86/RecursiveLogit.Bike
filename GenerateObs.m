global incidenceFull; 
    global Mfull;
    global Ufull;
    global Atts;
    global Op;
    global isLinkSizeInclusive;
    
    % Generate Obs
    % ---------------------------------------------------- 
    %x0=[-2.2825;-0.9675;-0.9392;-3.2907;-1.6103;0.7136;1.8116;0.9465;-5.3697;2.9506;-0.2985;0.1629;-0.3175;1.4226] %parameters
    x0=[-2.28116083;-0.97003141;-0.94020870;-3.29392602;-1.61022927;0.71264944;1.81025458;0.94731958;-5.36749695;2.94739261;-0.29814004;0.16450869;-0.31756774;1.42237100];
    loadData_2WaysNetwork;
    Op.n=length(x0);
    Op.x=x0;
    
    mu = 1; % MU IS NORMALIZED TO ONE
    % Parameter for the random term
    location = 0;
    scale = 1;
    euler = 0.577215665;  
    [lastIndexNetworkState,maxDest] = size(incidenceFull);
    [p q] = size(incidenceFull);
    nbobs = size(Obs,1);
    nbobsOD=ones(1,nbobs); %number of simulated obs to generate for each OD
      
    Mfull = getM(x0);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(x0);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    % Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    % Initialize
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    for i = 1:Op.n
        AttLc(i).value =  (Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end

    % b matrix:
    N = size(M,1);
    b = sparse(zeros(N,1));
    b(N) = 1;
    B = sparse(zeros(N, maxDest - lastIndexNetworkState));
    B(N,:) = ones(1,maxDest - lastIndexNetworkState);
    for i = 1: maxDest - lastIndexNetworkState
        B(1:lastIndexNetworkState,i) = Mfull(:, i+lastIndexNetworkState);
    end
    A = speye(size(M)) - M;
    Z = A\B;
    % Check feasible
    minele = min(Z(:));
    expVokBool = 1;
    if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
    end 
    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE? 
    D = (A * Z - B);
    resNorm = norm(D(:));
    if resNorm > OptimizeConstant.RESIDUAL
       expVokBool = 0;
    end
    if (expVokBool == 0)
            LL = OptimizeConstant.LL_ERROR_VALUE;
            grad = ones(Op.n,1);
            disp('The parameters not fesible')
    end
    
    IregularNetwork = incidenceFull(1:lastIndexNetworkState,1:lastIndexNetworkState);    
    dummy = lastIndexNetworkState + 1;
        
    % Obs = sparse(zeros(nbobs * nbobsOD, dummy));
    % Loop over all OD pairs
    timesampling=zeros(1,40);
    timeloop=zeros(1,40);
    for n = 1:nbobs
        n
        tic
        dest = Obs(n, 1);
        orig = Obs(n, 2);            
        addColumn = incidenceFull(:,dest);
        Incidence = IregularNetwork;
        Incidence(:,lastIndexNetworkState+1) = addColumn;
        expV = Z(:,dest - lastIndexNetworkState);
        expV = full(abs(expV));  
        V = log(expV);
        % Now we have all utilities, time to generate the observations
        for i = 1: nbobsOD
            SimObs((n-1)*nbobsOD + i, 1) = dest;
            SimObs((n-1)*nbobsOD + i, 2) = orig;
            k = orig;
            t = 3;
            while k ~= dummy
                ind = find(Incidence(k,:));
                nbInd = size(ind,2);
                bestUtilities = -1e6;
                for j = 1: nbInd
                    utility = U(k,ind(j)) + V(ind(j)) + random('ev',location,scale) - euler ;                   
                    if bestUtilities < utility
                        bestInd = ind(j);
                        bestUtilities = utility;
                    end
                end
                if bestInd ~= dummy
                    SimObs((n-1)*nbobsOD + i, t) = bestInd;
                    t = t + 1;
                end
                k = bestInd;
            end
            SimObs((n-1)*nbobsOD + i, t) = dest;
        end
        timeloop(n)=toc;
    end
