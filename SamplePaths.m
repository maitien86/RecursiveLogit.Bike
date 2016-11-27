global incidenceFull; 
global Mfull;
global LSatt;
global Atts;
global ODmat

loadDataForLinkFLows; 

%Parameters
nbOD=666;
%betaold:
%beta=[-1.94,-0.93,-0.85,-3.8,-1.6,0.33,1.49,0.7,-5.1,4.07,-0.22,0.34,-0.31,1.37];
%beta RL (old):
%beta=[-2.26805775 ,-1.00702583 ,-0.95312159, -3.29137429, -1.61254607, 0.71742112, 1.79588938 ,0.93605082,-5.42521752, 3.00972377, -0.30722496, 0.23232475,-0.32469659 ,1.42675390 ];
%beta RL (after left/right correction)
beta=[-2.25468977;-0.81323812;-1.01016929;-3.23739753;-1.60591263;0.73642977;1.80227092;0.92017333;-5.40745425;2.83432999;-0.28145550;-0.27620827;-1.83708946;1.36961231]; 
%betacommute:
%beta=[-3.03216802,-0.98362915 ,-1.30637436 ,-4.44948290 ,-1.52093996 ,0.89667749 ,1.76414608 ,1.02495418 ,-4.87537505 ,1.66597041 ,-0.28891665 ,0.29849662 ,-0.38794887 ,1.33599360 ];
%beta = [-5.2570;-4.6702;-5.0949;-1.5374;-5.1486]';
mu = 1; % MU IS NORMALIZED TO ONE
location = 0;
scale = 1;
euler = 0.577215665;  
r=20; %number of sampled path per OD pair

[nbStartLinks, nbAbsorbLinks] = size(incidenceFull);
lastIndexNetworkState = nbStartLinks;
u = beta(1) * Atts(1).value;
for i = 2:14 %change here for number of parameters
    u = u + beta(i) * Atts(i).value; 
end
expM = u;
expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
Mfull = expM;    
    
[p q] = size(Mfull);
M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState); 
M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));

% b matrix:
N = size(M,1);
b = sparse(zeros(N,1));
b(N) = 1;
B = sparse(zeros(N, nbAbsorbLinks - lastIndexNetworkState));
B(N,:) = ones(1,nbAbsorbLinks - lastIndexNetworkState);
for i = 1: nbAbsorbLinks - lastIndexNetworkState
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
    disp('The parameters not fesible')
end

IregularNetwork = incidenceFull(1:lastIndexNetworkState,1:lastIndexNetworkState);    
dummy = lastIndexNetworkState + 1;
[row,col]=find(ODmat);
Flow=zeros(nbAbsorbLinks,1);
start=tic;
for n = 1:length(row)
    n
    dest = col(n);
    orig = row(n);            
    addColumn = incidenceFull(:,dest);
    Incidence = IregularNetwork;
    Incidence(:,lastIndexNetworkState+1) = addColumn;
    expV = Z(:,dest - lastIndexNetworkState);
    expV = full(abs(expV));  
    V = log(expV);
    % Now we have all utilities, time to generate the observations
    nbDemand=ODmat(orig,dest)*r;
    for i = 1: nbDemand
        Obs((n-1)*nbDemand + i, 1) = dest;
        Obs((n-1)*nbDemand + i, 2) = orig;
        k = orig;
        t = 3;
        while k ~= dummy
            ind = find(Incidence(k,:));
            nbInd = size(ind,2);
            bestUtilities = -1e6;
            for j = 1: nbInd
                utility = u(k,ind(j)) + V(ind(j)) + random('ev',location,scale) - euler ;                                  
                if bestUtilities < utility
                    bestInd = ind(j);
                    bestUtilities = utility;
                end
            end
            if bestInd ~= dummy
                Obs((n-1)*nbDemand + i, t) = bestInd;
                Flow(bestInd)=Flow(bestInd)+1/r;
                t = t + 1;
            end
            k = bestInd;
        end
        Obs((n-1)*nbDemand + i, t) = dest;
    end
end
ElapsedTime=toc(start);
fprintf('\n Simulation time %d \n', ElapsedTime);

save('FlowSimMeth1r20.mat','Flow')
% save('Obs.mat','Obs')