% global incidenceFull; 
% global Mfull;
% global LSatt;
% global Atts;
% global ODmatCommute;
% global ODmatNonCommute;
 
loadDataPortland; %check here for restriction on links
ODmat=ODmatCommute+ODmatNonCommute;

% %Parameters
%betacommute
%beta=[-3.03216802,-0.98362915 ,-1.30637436 ,-4.44948290 ,-1.52093996 ,0.89667749 ,1.76414608 ,1.02495418 ,-4.87537505 ,1.66597041 ,-0.28891665 ,0.29849662 ,-0.38794887 ,1.33599360 ];
mu=1;
%beta modified (length-1)
%beta=[-3.26805775 ,-1.00702583 ,-0.95312159, -3.29137429, -1.61254607, 0.71742112, 1.79588938 ,0.93605082,-5.42521752, 3.00972377, -0.30722496, 0.23232475,-0.32469659 ,1.42675390 ];
%beta invented
beta=[-4.0 ,-1.00702583 ,-0.95312159, -3.29137429, -1.61254607, 0.5, 0.8 ,0.7,-5.42521752, 3.00972377, -0.30722496, 0.23232475,-0.32469659 ,1.42675390 ];


[nbStartLinks, nbAbsorbLinks] = size(incidenceFull);
nbOD=nbAbsorbLinks-nbStartLinks;
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
%additionned flows to all destination
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
        [expV, expVokBool] = getExpV(M); % vector with value functions for given beta 
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
    
%     %% Add correction term               
%     ips(size(incidenceFull,2),1) = 0;      
%     I = find(incidenceFull);
%     [nbnonzero, c] = size(I);
%     ind1 = zeros(nbnonzero,1);
%     ind2 = zeros(nbnonzero,1);
%     s = zeros(nbnonzero,1);
%     loop=tic;
%     for i = 1:nbnonzero
%         [k a] = ind2sub(size(incidenceFull), I(i));
%         ind1(i) = k;
%         ind2(i) = a;
%         s(i) =  ips(a);
%     end 
%     timeloop(n)=toc(loop);
%     createlsatt=tic;
%     lsatt = sparse(ind1, ind2, s, size(incidenceFull,1), size(incidenceFull,2));        
%     LSatt(n).value = lsatt;
%     timecreatelsatt(n)=toc(createlsatt);
%     timeremaining(n)=toc(remaining);
end
SimTime=toc(startSimTime);

%% Create table
TableFlows=table(sumF(1:308520),linkNumber,'VariableNames',{'LinkIDs','SimulatedFlow'});
%writetable(TableFlows,'LinkFlowsPortland.csv');

%% Scatter Plot
SimFlow=sumF(1:308520);
VolTot=linkAttributes(:,27);

%density scatter plot
values=hist3([VolTot,SimFlow],[200,200]); %change for more resolution
%values=values(1:200,1:200); %to zoom in on lower corner
figure
imagesc(values,[0,20])
colorbar
colormap jet
axis xy
xticklabels = 0:820:4100;
xticks = linspace(1, 200, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:700:2800;
yticks = linspace(1, 200, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
xlabel('RL model link flows')
ylabel('Path-based model link flows')

%scatter plot
figure
scatter(SimFlow,VolTot);
x=[1:3000];
hold on
plot(x,x,'k')
hold on
rectangle('Position',[0 0 1000 750],'LineWidth',3)
xlabel('RL model link flows')
ylabel('Path-based model link flows')

%% Analyzing bike facilities
indicesBikeLane=linkAttributes(:,19)==1;
indicesBikeBlvd=linkAttributes(:,15)==1;
indicesRMUP=ismember(FacType,0);

figure
scatter(SimFlow(indicesBikeLane),VolTot(indicesBikeLane),'MarkerFaceColor','b','MarkerEdgeColor','k')
legend('Bike Lanes','Location','southeast')
x=[1:3000];
hold on
plot(x,x,'k');
xlabel('RL model link flows')
ylabel('Path-based model link flows')
figure
scatter(SimFlow(indicesBikeBlvd),VolTot(indicesBikeBlvd),'MarkerFaceColor','g','MarkerEdgeColor','k')
legend('Bike Boulevards','Location','southeast')
x=[1:2000];
hold on
plot(x,x,'k');
xlabel('RL model link flows')
ylabel('Path-based model link flows')
figure
scatter(SimFlow(indicesRMUP),VolTot(indicesRMUP),'MarkerFaceColor','r','MarkerEdgeColor','k')
legend('RMUPs','Location','southeast')
x=[1:3000];
hold on
plot(x,x,'k');
xlabel('RL model link flows')
ylabel('Path-based model link flows')

%% Compute r squared
x=CompareFlowsSorted(:,2);%VolTot(indicesBikeLane);
y=CompareFlowsSorted(:,3);%SimFlow(indicesBikeLane);
average=sum(x)/length(x);
averagevect=average*ones(length(x),1);
DiffVolTot=x-averagevect;
SSTot=sum(DiffVolTot.*DiffVolTot);
Diff=x-y;
SSRes=sum(Diff.*Diff);
RSquared=1-(SSRes/SSTot);