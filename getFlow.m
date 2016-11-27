%   Compute Link Size attribute from data
%   
%%
function [ok, Flow] = getFlow(beta)
    global incidenceFull; 
    global Obs;     % Observation
    global Atts;
    %global EstimatedTime;
    %global TurnAngles;
    %global LeftTurn;
    %global Uturn;    
    % ----------------------------------------------------
    mu = 1; % MU IS NORMALIZED TO ONE
    nParamsFlow = length(beta);
    Obs_temp = Obs;
    Obs_temp(find(Obs)) =  Obs(find(Obs)) + 1;
    temp_Atts = objArray(nParamsFlow);
    % Save data
    for i = 1:nParamsFlow
        temp_Atts(i).value = Atts(i).value;
    end
    IF = (incidenceFull);    
    incidenceFull = transferMarix(incidenceFull);
    for i = 1:nParamsFlow
        temp_Atts(i).value = sparse(transferMarix(Atts(i).value));
    end
    lastIndexNetworkState = size(incidenceFull,1);
    dest = unique(Obs_temp(:,1));
    orig = unique(Obs_temp(:,2));
    incidenceFull(1,orig') = 1;
    temp_Atts(1).value(1,orig') = 0.05;
    incidenceFull(dest,lastIndexNetworkState) = 1;
    temp_Atts(1).value(dest,lastIndexNetworkState) = 0.05;    
    %M = getStandardM(beta,false);
    u = beta(1) * temp_Atts(1).value;
    for i = 2:nParamsFlow
        u = beta(i) * temp_Atts(i).value + u;
    end
    u = sparse(u);
    expM = u ;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    M = incidenceFull .* expM;
    [expV, expVokBool] = getExpV(M);
    if (expVokBool == 0)
       ok = false;
       Flow = [];
       disp('ExpV is not fesible')
       return; 
    end  
    P = getP(expV, M);
    G = sparse(zeros(size(expV)));
    G(1) = 1;
    I = speye(size(P));
    F = (I-P')\G;                        
    if (min(F) < 0)                
        ToZero = find(F <= 0);
        for i=1:size(ToZero,1)
            F(ToZero(i)) = 1e-9;
        end
    end
    ips = F;
    ips = ips(2:lastIndexNetworkState - 1 ); 
    Flow = ips;
    % Recover real size
    incidenceFull = IF;
    ok = true;
end