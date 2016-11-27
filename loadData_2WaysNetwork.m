%%   Load Bike Route Choice Data
% Loading data for the 2 ways network
%%
disp('Loading data ...')

global file_linkIncidence;
global file_observations;
global file_linkAttributes;
global incidenceFull;
global linkAttributes;
global Obs;
global nbobs; 
global Atts;
global AttsNames;
global LSatt;

file_linkIncidence = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/Incidence_2WaysNetwork.txt'; %matrix of all possible state transitions (1 if transition between 2 states is possible, 0 else)
file_observations = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/TripsObservations_2WaysNetwork.txt';
file_linkAttributes = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/fullLinksAtts_2WaysNetwork.txt';
file_linkAttributes2 = '/home/zimmae/Documents/RL2.6.DeC.GLS/Input/LinkIDs.xlsx';

incidenceFull = spconvert(load(file_linkIncidence));
linkAttributes = load(file_linkAttributes) ; %each link results in two states! (two ways)
Obs = spconvert(load(file_observations));
TableAttributes = readtable(file_linkAttributes2);

%temporary here
%usually no need to remove obs, they are all valid since links are 2 ways
% FalseObs=FindFalseObs();
% [indexFalseObsRow,indexFalseObsCol]=find(FalseObs);
% indexFalseObs=unique(indexFalseObsRow);
% Obs(indexFalseObs,:)=[];

[nbobs, ~] = size(Obs);
[nbRealStates, nbTotalStates] = size(incidenceFull); %The states are all possible links
fprintf('taille de incidence %d et %d \n',nbRealStates,nbTotalStates);
linkLength =  linkAttributes(:,6);
linkNumber = linkAttributes(:,1);
FromNode=linkAttributes(:,2);
ToNode=linkAttributes(:,3);

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
AngleFromNode(nbTotalStates)=0;
AngleToNode(nbTotalStates)=0;
NoTurnDummy = sparse(nbRealStates,nbTotalStates);
CrossingLeft = sparse(nbRealStates,nbTotalStates);
GoingRight = sparse(nbRealStates,nbTotalStates);
GoingRightPossible = sparse(nbRealStates,nbTotalStates);
Intersection = sparse(nbRealStates,nbTotalStates);
[row, col]=find(incidenceFull>0);
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
bikeFacility=table2cell(TableAttributes(:,17));%bike facility type
bothways=table2cell(TableAttributes(:,5));%dummy for whether the link can be taken both ways
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
stop = table2cell(TableAttributes(:,29));
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
signal = table2cell(TableAttributes(:,30));
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

%Speed limit
speedlimit =linkAttributes(:,38);
NaN=find(isnan(speedlimit));
speedlimit(NaN)=0; %speed limit (between 5 and 65 mph)
speedlimitDummy=speedlimit>50; %dummy for speed limit higher than 50mph (=80kmh)

linkLength(nbTotalStates) = 0;
linkDzp(nbTotalStates) = 0;
linkDzn(nbTotalStates) = 0;
diffDz(nbTotalStates) = 0;
slope(nbTotalStates) = 0;
upslope(nbTotalStates) = 0;
slopeDummy6p(nbTotalStates)=0;
upslopeDummy4p(nbTotalStates)=0;
upslopeDummy4p6p(nbTotalStates)=0;
estVolumeLargeDummy8000(nbTotalStates)=0;
estVolumeSmallDummy(nbTotalStates)=0;
estVolumeHugeDummy20000(nbTotalStates)=0;
BikeFacilityDummy(nbTotalStates)=0;
stopdummy(nbTotalStates)=0;
signaldummy(nbTotalStates)=0;
nosignaldummy(nbTotalStates)=0;
speedlimitDummy(nbTotalStates)=0;
BikeBlvd(nbTotalStates)=0;
Rmup(nbTotalStates)=0;
BikeLane(nbTotalStates)=0;
BridgeNoFac(nbTotalStates)=0;
BridgeWithFac(nbTotalStates)=0;
Bridge(nbTotalStates)=0;

icdLength = incidenceFull * spdiags(linkLength,0,nbTotalStates,nbTotalStates);
icdDzp = incidenceFull * spdiags(linkDzp,0,nbTotalStates,nbTotalStates);
icdDzn = incidenceFull * spdiags(linkDzn,0,nbTotalStates,nbTotalStates);
icdDiffDz =  incidenceFull * spdiags(diffDz,0,nbTotalStates,nbTotalStates);
icdupslopeDummy4p = incidenceFull * spdiags(upslopeDummy4p,0,nbTotalStates,nbTotalStates);
icdupslopeDummy4p6p = incidenceFull * spdiags(upslopeDummy4p6p,0,nbTotalStates,nbTotalStates);
icdestVolumeLargeDummy8000 = incidenceFull * spdiags(estVolumeLargeDummy8000,0,nbTotalStates,nbTotalStates);
estVolumeLargeStateDummy = spdiags(estVolumeLargeDummy8000,0,nbRealStates,nbRealStates) * incidenceFull;
estVolumeHugeStateDummy = spdiags(estVolumeHugeDummy20000,0,nbRealStates,nbRealStates) * incidenceFull;
icdestVolumeSmallDummy = incidenceFull * spdiags(estVolumeSmallDummy,0,nbTotalStates,nbTotalStates);
icdestVolumeHugeDummy20000 = incidenceFull * spdiags(estVolumeHugeDummy20000,0,nbTotalStates,nbTotalStates);
icdBikeFacilityDummy = incidenceFull * spdiags(BikeFacilityDummy,0,nbTotalStates,nbTotalStates);
icdStopDummy = incidenceFull * spdiags(stopdummy,0,nbTotalStates,nbTotalStates);
icdSignalDummy = incidenceFull * spdiags(signaldummy,0,nbTotalStates,nbTotalStates);
SignalStateDummy = spdiags(signaldummy,0,nbRealStates,nbRealStates) *incidenceFull;
noSignalStateDummy = spdiags(nosignaldummy,0,nbRealStates,nbRealStates) *incidenceFull;
icdspeedlimitDummy = incidenceFull * spdiags(speedlimitDummy,0,nbTotalStates,nbTotalStates);
icdBikeBlvd = incidenceFull * spdiags(BikeBlvd,0,nbTotalStates,nbTotalStates);
icdRmup = incidenceFull * spdiags(Rmup,0,nbTotalStates,nbTotalStates);
icdBikeLane = incidenceFull * spdiags(BikeLane,0,nbTotalStates,nbTotalStates);
icdBridgeNoFac = incidenceFull * spdiags(BridgeNoFac,0,nbTotalStates,nbTotalStates);
icdBridgeWithFac = incidenceFull * spdiags(BridgeWithFac,0,nbTotalStates,nbTotalStates);
icdBridge = incidenceFull * spdiags(Bridge,0,nbTotalStates,nbTotalStates);
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
%load('/home/zimmae/Documents/RL2.6.DeC.GLS/LSatt2WaysNetwork.mat'); %next times load simply LSatt

% %get global link size:
% beta = [-5.2570;-4.6702;-5.0949;-1.5374;-5.1486]';
% [ok,LS]=getFlow(beta);
% LS(nbTotalStates)=0;
% icdLS=incidenceFull*spdiags(LS,0,nbTotalStates,nbTotalStates);
% Atts(15).value = icdLS;
