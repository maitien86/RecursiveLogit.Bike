global incidenceFull; 
global Mfull;
global LSatt;
global Atts;
global ODmat

%Load network attributes and OD demand matrix
loadDataForLinkFLows;

%Parameters
nbOD=666;

%Choice of beta values
%beta RL: 
beta=[-2.25468977;-0.81323812;-1.01016929;-3.23739753;-1.60591263;0.73642977;1.80227092;0.92017333;-5.40745425;2.83432999;-0.28145550;-0.27620827;-1.83708946;1.36961231]; 
%beta RL-LS:
%beta=[-2.28388074;-0.82069841;-1.01628239;-3.15143101;-1.59705910;0.76107886;1.81426372;0.86818907;-4.55540630;1.99111811;-0.29179754;-0.33991884;-1.86365984;1.33100569;-0.23549472];
%betacommute:
%beta=[-3.03216802,-0.98362915 ,-1.30637436 ,-4.44948290 ,-1.52093996 ,0.89667749 ,1.76414608 ,1.02495418 ,-4.87537505 ,1.66597041 ,-0.28891665 ,0.29849662 ,-0.38794887 ,1.33599360 ];

mu = 1; % MU IS NORMALIZED TO ONE

%% Link flow computation
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

timesystemlin=zeros(1,nbOD);
timegetV=zeros(1,nbOD);
timegetP=zeros(1,nbOD);
timegetPold=zeros(1,nbOD);
timegetF=zeros(1,nbOD);

sumF=zeros(lastIndexNetworkState,1);
ProblemWithFlow=0;
F=zeros(lastIndexNetworkState+1,nbOD);

startSimTime=tic;
for n = 1:nbOD 
    n
    start=tic;
    destlink = lastIndexNetworkState+n;
    if true %expV0(dest) == 0   
        M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,destlink);
        Vtime=tic;
        [expV, expVokBool] = getExpVLinkFlows(M); % vector with value functions for given beta 
        timegetV(n)=toc(Vtime);
        if (expVokBool == 0)
            ExpV_is_ok = false;
            disp('ExpV is not feasible')
            return; 
        end  
        %P = getP_old(expV, M);
        Ptime=tic;
        P=getP(expV,M);
        timegetP(n)=toc(Ptime);
    end
    %G = sparse(zeros(size(expV)));
    %G(orig) = 1;
    G = ODmat(:,destlink);
    G(lastIndexNetworkState+1)=0;
    I = speye(size(P));
    system=tic;
    F(:,n) = (I-P')\G;   
    timesystemlin(n)=toc(system);
    timegetF(n)=toc(start);
    if (min(F(:,n)) < 0)                
        ToZero = find(F(:,n) <= 0);
        for i=1:size(ToZero,1)
            F(ToZero(i),n) = 1e-9;
        end
    end
%    remaining=tic;
    ips = F(:,n);
    ips = ips(1:lastIndexNetworkState); % dummy link should not be included
    if any(isnan(ips))
        ProblemWithFlow=ProblemWithFlow+1;
    else
        sumF = sumF+ips;
    end
end
SimTime=toc(startSimTime);
fprintf('\n Simulation time %d \n', SimTime);

%Uncomment to save vector of link flows
%save('FlowSystem.mat','sumF')

%% Assign simulated flow to Link ID
LinkFlows=zeros(nbRealStates,2);
uniqueLinkNumber=unique(linkNumber);
LinkFlowsRealNetwork=zeros(length(uniqueLinkNumber),2);
for i = 1:nbRealStates
    LinkFlows(i,1)=linkNumber(i);
    LinkFlows(i,2)=sumF(i);
end
%Merge flows of states corresponding to same link (Way to and from)
for i = 1:length(uniqueLinkNumber)
    LinkFlowsRealNetwork(i,1)=linkNumber(i);
    LinkFlowsRealNetwork(i,2)=sum(LinkFlows(linkNumber==uniqueLinkNumber(i),2));
end

%% Save to file for QGIS
TableFlows=table(LinkFlowsRealNetwork(:,1),LinkFlowsRealNetwork(:,2),'VariableNames',{'LinkIDs','SimulatedFlow'});
TableFlows=[TableFlows TableAttributes(:,9)];

%uncomment to save the table
%writetable(TableFlows,'FlowAllLinksWithBikeRest.csv');


% COMPARISON WITH OBSERVED LINK COUNTS----------------------------------
%% Read link counts
file_linkCounts = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/Daily_for_Emma.xlsx';
TableCounts = readtable(file_linkCounts);
LinkIDcount=table2cell(TableCounts(:,25));
DailyCount=table2array(TableCounts(:,4));
YearCount=table2array(TableCounts(:,7));
Del=ismember(LinkIDcount,'NA');
LinkIDcount(Del)=[];
LinkIDcount=str2double(LinkIDcount);
DailyCount(Del)=[];
YearCount(Del)=[];
uniqueLinkID=unique(LinkIDcount);
Counts=zeros(length(uniqueLinkID),2);
for i = 1:length(uniqueLinkID)
    logicalindex=LinkIDcount==uniqueLinkID(i);
    Counts(i,1)=uniqueLinkID(i);
    Counts(i,2)=sum(DailyCount(logicalindex))/sum(logicalindex);
end

%% Read flows from other models
file_allmodels='/home/zimmae/Documents/RL2.6.DeC.GLS/Input/AllModels_AvgDaily_nb.xlsx';
LinkCountAllModels=readtable(file_allmodels,'range','A2:B18842');
LinkCountAllModels=table2array(LinkCountAllModels);

%% Comparison with counts
CompareFlows=zeros(length(uniqueLinkID),6);
CompareFlows(:,1:2)=Counts;
for i = 1:length(uniqueLinkID)
    CompareFlows(i,3)=sum(LinkFlows(linkNumber==uniqueLinkID(i),2));
    if isempty(find(LinkCountAllModels(:,1)==uniqueLinkID(i)))
        CompareFlows(i,4)=0;
    else
        CompareFlows(i,4)=LinkCountAllModels(LinkCountAllModels(:,1)==uniqueLinkID(i),2);
    end
end
CompareFlows(:,3)=CompareFlows(:,3); %normalize so it sums to observed counts
for i = 1:length(uniqueLinkID)
    CompareFlows(i,5)=abs(CompareFlows(i,2)-CompareFlows(i,3));
    CompareFlows(i,6)=abs(CompareFlows(i,3)-CompareFlows(i,4));
    CompareFlows(i,7)=abs(CompareFlows(i,2)-CompareFlows(i,4));
end
[Y,I]=sort(CompareFlows(:,5));
CompareFlowsSorted=CompareFlows(I,:);

%% Prepare xcel table for Emma
TheIDs=CompareFlowsSorted(:,1);
indices=zeros(length(TheIDs),1);
for i=1:length(TheIDs)
    indices(i)=find(linkNumber(1:21192)==TheIDs(i));
end

TableOutput=table(TheIDs,'VariableNames',{'LinkIDs'});
TableOutput2=table(CompareFlowsSorted(:,2),CompareFlowsSorted(:,3),CompareFlowsSorted(:,5),'VariableNames',{'SimulatedFlow','ObservedCount','Difference'});
TableOutput=[TableOutput TableAttributes(indices,9) TableAttributes(indices,11) TableOutput2];

%% Scatter Plot
untest=bikeFacility(indices);
rmupindices=strcmp(untest,'rmup');
bikelaneindices1=strcmp(untest,'bike lanes');
bikelaneindices2=strcmp(untest,'bike lane FT');
bikelaneindices3=strcmp(untest,'bike lane TF');
bikelaneindices=bikelaneindices1 | bikelaneindices2 | bikelaneindices3;
otherindices=~(rmupindices | bikelaneindices);
figure
scatter(CompareFlowsSorted(rmupindices,2),CompareFlowsSorted(rmupindices,3),'MarkerFaceColor','r');
hold on
scatter(CompareFlowsSorted(bikelaneindices,2),CompareFlowsSorted(bikelaneindices,3),'MarkerFaceColor','b');
hold on
scatter(CompareFlowsSorted(otherindices,2),CompareFlowsSorted(otherindices,3),'MarkerFaceColor','k');
x=[1:2000];
hold on
plot(x,x);
xlabel('Observed link count')
ylabel('Simulated link flow')
legend('RMUPs','Bike lanes','other links')