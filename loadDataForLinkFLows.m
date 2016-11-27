%%   Load Bike Route Choice Data
% Loading data for the 2 ways network
%%
disp('Loading data ...')

global file_linkIncidence;
global file_linkAttributes;
global incidenceFull;
global linkAttributes;
global Atts;
global AttsNames;
global ODmat;
global LSatt;

file_linkIncidence = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/IncidenceForLinkFLows.txt';
file_linkAttributes = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/fullLinksAttsForLinkFlows.txt';
file_linkAttributes2 = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/BikeFacility.xlsx';
file_ODMat = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/ODmat.txt';

%load incidence matrix
incidenceFull = spconvert(load(file_linkIncidence));
%load network attributes
linkAttributes = load(file_linkAttributes) ; %each link results in two states! (two ways)
%load OD demand matrix
ODmat=spconvert(load(file_ODMat));
%load additional network attributes
TableAttributes = readtable(file_linkAttributes2);

%define number of states
[nbStatesRealStart,nbStatesRealStartAbsorb] = size(incidenceFull);
nbRealStates = length(linkAttributes(:,1));
fprintf('taille de incidence %d et %d \n',nbStatesRealStart,nbStatesRealStartAbsorb);

%Basic link attributes
linkLength =  linkAttributes(:,6);
linkNumber = linkAttributes(:,1);
FromNode=linkAttributes(:,2);
ToNode=linkAttributes(:,3);
BikeRest=linkAttributes(:,14);

%Implement bike restriction
%Too high length for links with bike restriction
BikeRestInd=find(BikeRest);
linkLength(BikeRestInd)=20090; %maximum length
linkLength(10898)=80000;
linkLength(32029)=80000;
linkLength(10837)=80000;
linkLength(32090)=80000;

%Elevation
linkDzp = linkAttributes(:,7); %elevation rise
linkDzn = linkAttributes(:,8); %elevation fail
diffDz = linkDzp - linkDzn; %rise-fail

%Slope
slope = diffDz./linkLength; %average slope
upslope = linkDzp./linkLength; %average upslope
slopeDummy6p = (slope>=0.06); %dummy for average slope higher than 6%
upslopeDummy4p = (upslope>=0.04); %dummy for average upslope higher than 4%
upslopeDummy4p6p = (upslope>=0.04) & (upslope < 0.06); %dummy for average upslope between 4 and 6%

%Bridge
%BridgeWF = linkAttributes(:,27); %seems wrong in the data
%BridgeNF = linkAttributes(:,28); %seems wrong in the data
Bridge = linkAttributes(:,33);

%Turn attributes 
AngleFromNode = linkAttributes(:,10);
AngleToNode = linkAttributes(:,11);
AngleFromNode(nbStatesRealStartAbsorb)=0;
AngleToNode(nbStatesRealStartAbsorb)=0;
NoTurnDummy = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
CrossingLeft = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
GoingRight = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
GoingRightPossible = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
Intersection = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
[row, col]=find(incidenceFull(nbRealStates,nbRealStates)>0);
L=length(row);
for k=1:L
    NoTurnDummy(row(k),col(k))=goingStraight(AngleToNode(row(k)),AngleFromNode(col(k)));
end
for k=1:L
    CrossingLeft(row(k),col(k))=goingLeft(AngleToNode(row(k)),AngleFromNode(col(k)));
    GoingRight(row(k),col(k)) = goingRight(AngleToNode(row(k)),AngleFromNode(col(k)));
end
for k=1:nbRealStates
    if any(GoingRight(k,:))%if link k has any outgoing link going right
        GoingRightPossible(k,:)=incidenceFull(k,:); %link k has the possibility to go right 
        if any(NoTurnDummy(k,:))% if link k has any outgoing link going straight
            if any(CrossingLeft(k,:))
                Intersection(k,:)=incidenceFull(k,:);%link k has also the possibility to go straight and left so link k is at an intersection
            end
        end
    end
end

CrossingLeft = CrossingLeft.*Intersection; %We only consider left turns at intersections where there is also an option to go right and straight
GoingRight = GoingRight.*Intersection;
icdNoTurnIntersection=NoTurnDummy.*Intersection;
icdNoTurn=NoTurnDummy;

%Volume
estVolume = linkAttributes(:,26);%estimated traffic volume
NaN = find(isnan(estVolume));
estVolume(NaN)=0;
estVolumeLargeDummy8000=estVolume>=8000 & estVolume<20000; %dummy for traffic volume greater than 8000 = proxy for highway or very large road
estVolumeSmallDummy = estVolume < 1000; %dummy for small traffic volume
estVolumeHugeDummy20000=estVolume>=20000; %dummy for traffic volume greater than 20000

%Bike Facility & One Way restrictions
bikeFacility=table2cell(TableAttributes(:,11));%bike facility type
bothways=table2cell(TableAttributes(:,13));%dummy for whether the link can be taken both ways
bothwaysdummy=ismember(bothways,'');%dummy for whether the link can be taken both ways
zeroindex=find(bothwaysdummy==0);
FTrestriction=ismember(bothways,'FT');%dummy for whether there is a restriction in the FT direction

%Specific bike facilities
BikeBlvdWay1=ismember(bikeFacility,'blvd');
BikeBlvdWay2=BikeBlvdWay1;
BikeBlvd=[BikeBlvdWay1;BikeBlvdWay2];%dummy for Bike Boulevard (=street only for bikes)
RmupWay1=ismember(bikeFacility,'rmup');
RmupWay2=RmupWay1;
Rmup=[RmupWay1;RmupWay2];%dummy for regional multi use path (=small street for pedestrian and bikers)
BikeLaneBothWays=ismember(bikeFacility,'bike lanes');
BikeLaneWay1=ismember(bikeFacility,'bike lane FT');
BikeLaneWay2=ismember(bikeFacility,'bike lane TF');
BikeLaneWay1=BikeLaneWay1 + BikeLaneBothWays;
BikeLaneWay2=BikeLaneWay2 + BikeLaneBothWays;
%account for the fact that "Way1" includes some TF links too (the ones with
%FT restriction)
intersect1=BikeLaneWay1.*FTrestriction;%there shouldn't be any link that has a FT restriction but also a Bike Lane in the FT dir.
intersect2=BikeLaneWay2.*FTrestriction;%but here there could be links with FT restriction and Bike Lane in TF direction
BikeLaneWay1=BikeLaneWay1-intersect1+intersect2;
BikeLaneWay2=BikeLaneWay2-intersect2+intersect1;
BikeLane=[BikeLaneWay1;BikeLaneWay2];%dummy for Bike Lane (=lane for bikers on the side of the road)

%Bridges with/without bike facilities
BikeFacilityDummy = BikeBlvd | Rmup | BikeLane;
noBikeFacilityDummy = ones(length(BikeFacilityDummy),1)-BikeFacilityDummy;
BridgeNoFac = Bridge.*noBikeFacilityDummy;
BridgeWithFac = Bridge.*BikeFacilityDummy;

%Stops
stop = table2cell(TableAttributes(:,24));
stopWay1= ismember(stop,'FT');
stopWay2=ismember(stop,'TF');
stopbothways=ismember(stop,'Both');
stopWay1 = stopWay1 + stopbothways;
stopWay2 = stopWay2 + stopbothways;
%account for the fact that "Way1" includes some TF links too (the ones with
%FT restriction)
intersect1=stopWay1.*FTrestriction;
intersect2=stopWay2.*FTrestriction;
stopWay1=stopWay1-intersect1+intersect2;
stopWay2=stopWay2-intersect2+intersect1;
stopdummy = [stopWay1;stopWay2];%dummy for whether there is a stop at the end of the road

%signals
signal = table2cell(TableAttributes(:,25));
signalWay1=ismember(signal,'FT');
signalWay2=ismember(signal,'TF');
signalbothways=ismember(signal,'Both');
signalWay1=signalWay1 + signalbothways;
signalWay2=signalWay2 + signalbothways;
%account for the fact that "Way1" includes some TF links too (the ones with
%FT restriction)
intersect1=signalWay1.*FTrestriction;
intersect2=signalWay2.*FTrestriction;
signalWay1=signalWay1-intersect1+intersect2;
signalWay2=signalWay2-intersect2+intersect1;
signalWay2(zeroindex)=[];%removes links that don't have a "Way2" because they have a restriction in FT or TF direction
signaldummy=[signalWay1;signalWay2];%dummy for whether there is a signal at the end of the road
nosignaldummy=ones(length(signaldummy),1)-signaldummy;

linkLength(nbStatesRealStartAbsorb) = 0;
linkDzp(nbStatesRealStartAbsorb) = 0;
linkDzn(nbStatesRealStartAbsorb) = 0;
diffDz(nbStatesRealStartAbsorb) = 0;
slope(nbStatesRealStartAbsorb) = 0;
upslope(nbStatesRealStartAbsorb) = 0;
slopeDummy6p(nbStatesRealStartAbsorb)=0;
upslopeDummy4p(nbStatesRealStartAbsorb)=0;
upslopeDummy4p6p(nbStatesRealStartAbsorb)=0;
estVolumeLargeDummy8000(nbStatesRealStartAbsorb)=0;
estVolumeSmallDummy(nbStatesRealStartAbsorb)=0;
estVolumeHugeDummy20000(nbStatesRealStartAbsorb)=0;
BikeFacilityDummy(nbStatesRealStartAbsorb)=0;
stopdummy(nbStatesRealStartAbsorb)=0;
signaldummy(nbStatesRealStartAbsorb)=0;
nosignaldummy(nbStatesRealStartAbsorb)=0;;
BikeBlvd(nbStatesRealStartAbsorb)=0;
Rmup(nbStatesRealStartAbsorb)=0;
BikeLane(nbStatesRealStartAbsorb)=0;
BridgeNoFac(nbStatesRealStartAbsorb)=0;
BridgeWithFac(nbStatesRealStartAbsorb)=0;
Bridge(nbStatesRealStartAbsorb)=0;

icdLength = incidenceFull * spdiags(linkLength,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdDzp = incidenceFull * spdiags(linkDzp,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdDzn = incidenceFull * spdiags(linkDzn,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdDiffDz =  incidenceFull * spdiags(diffDz,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdupslopeDummy4p = incidenceFull * spdiags(upslopeDummy4p,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdupslopeDummy4p6p = incidenceFull * spdiags(upslopeDummy4p6p,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdestVolumeLargeDummy8000 = incidenceFull * spdiags(estVolumeLargeDummy8000,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
estVolumeLargeStateDummy = spdiags(estVolumeLargeDummy8000,0,nbStatesRealStart,nbStatesRealStart) * incidenceFull;
estVolumeHugeStateDummy = spdiags(estVolumeHugeDummy20000,0,nbStatesRealStart,nbStatesRealStart) * incidenceFull;
icdestVolumeSmallDummy = incidenceFull * spdiags(estVolumeSmallDummy,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdestVolumeHugeDummy20000 = incidenceFull * spdiags(estVolumeHugeDummy20000,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBikeFacilityDummy = incidenceFull * spdiags(BikeFacilityDummy,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdStopDummy = incidenceFull * spdiags(stopdummy,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdSignalDummy = incidenceFull * spdiags(signaldummy,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
SignalStateDummy = spdiags(signaldummy,0,nbStatesRealStart,nbStatesRealStart) *incidenceFull;
noSignalStateDummy = spdiags(nosignaldummy,0,nbStatesRealStart,nbStatesRealStart) *incidenceFull;
icdBikeBlvd = incidenceFull * spdiags(BikeBlvd,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdRmup = incidenceFull * spdiags(Rmup,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBikeLane = incidenceFull * spdiags(BikeLane,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBridgeNoFac = incidenceFull * spdiags(BridgeNoFac,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBridgeWithFac = incidenceFull * spdiags(BridgeWithFac,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBridge = incidenceFull * spdiags(Bridge,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
LargeTraffic = icdestVolumeLargeDummy8000 | estVolumeLargeStateDummy;
HugeTraffic = icdestVolumeHugeDummy20000 | estVolumeHugeStateDummy; %Large traffic on either link k or link a when going from k to a
CrossingLeftWithLargeTrafficUnsign = CrossingLeft.*LargeTraffic.*noSignalStateDummy;
CrossingLeftWithHugeTrafficUnsign = CrossingLeft.*HugeTraffic.*noSignalStateDummy;
CrossingLeftWithSignal = CrossingLeft.*SignalStateDummy;
icdNoTurnWithSignal = icdNoTurn.*SignalStateDummy;

Atts  = objArray(5);
Atts(1).value = icdLength ./ 1000;
Atts(2).value = (icdestVolumeLargeDummy8000.*icdLength) ./1000;
Atts(3).value = (icdestVolumeHugeDummy20000.*icdLength) ./1000;
%Atts(4).value = icdDzp ./ 100;
% Atts(2).value = (icdDzp .* Atts(1).value) ./ 100;
%Atts(4).value = (icdupslopeDummy4p6p.*icdLength) ./1000;
Atts(4).value = (icdupslopeDummy4p.*icdLength) ./1000;
Atts(5).value = incidenceFull;
Atts(6).value = (icdBikeBlvd.*icdLength) ./1000;
Atts(7).value = (icdRmup.*icdLength) ./1000;
Atts(8).value = (icdBikeLane.*icdLength) ./1000;
%Atts(8).value = (icdspeedlimitDummy.*icdLength)./1000;
%Atts(10).value = icdStopDummy;
%Atts(10).value = icdsignaldummy;
Atts(9).value = icdBridge;
Atts(10).value = icdBridgeWithFac;
Atts(11).value = icdNoTurnIntersection;
%Atts(13).value = CrossingLeft;
Atts(12).value = CrossingLeftWithLargeTrafficUnsign;
Atts(13).value = CrossingLeftWithHugeTrafficUnsign;
Atts(14).value = icdNoTurn;
%Atts(15).value = CrossingLeftWithSignal;

AttsNames = {varname(icdLength);varname(icdestVolumeLargeDummy8000);varname(icdestVolumeHugeDummy20000);varname(icdupslopeDummy4p);varname(incidenceFull);varname(icdBikeBlvd);varname(icdRmup);varname(icdBikeLane);varname(icdBridge);varname(icdBridgeWithFac);varname(icdNoTurnIntersection);varname(CrossingLeftWithLargeTrafficUnsign);varname(CrossingLeftWithHugeTrafficUnsign);varname(icdNoTurn)};

%get link Size attribute
%getLinkSizeAtt_new(); %The first time compute the linksize attribute with
%this function and save it
%load('/home/zimmae/Documents/RL2.6.DeC.GLS/LSatt.mat'); %next times load simply LSatt

% %get global link size:
% beta = [-5.2570;-4.6702;-5.0949;-1.5374;-5.1486]';
% [ok,LS]=getFlow(beta);
% LS(nbStatesRealAbsorb)=0;
% icdLS=incidenceFull*spdiags(LS,0,nbStatesRealAbsorb,nbStatesRealAbsorb);
% Atts(15).value = icdLS;