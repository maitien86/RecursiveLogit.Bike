%%   Load Bike Route Choice Data
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
global isLinkSizeInclusive;
global ODmatCommute;
global ODmatNonCommute;

file_linkIncidence = './Input/IncidencePortlandNetwork.txt'; %matrix of all possible state transitions (1 if transition between 2 states is possible, 0 else)
file_linkAttributes = './Input/fullLinksAttsPortlandNetwork.txt';
file_linkAttributes2 = './Input/LinkAttributesPortland.xlsx';
file_turnattributes = './Input/PortlandMetroBikeNetwork2010_turns.txt'
file_ODmatCommute = './Input/ODmatCommute.txt';
file_ODmatNonCommute = './Input/ODmatNonCommute.txt';

incidenceFull = spconvert(load(file_linkIncidence));
ODmatCommute=spconvert(load(file_ODmatCommute));
ODmatNonCommute=spconvert(load(file_ODmatNonCommute));
linkAttributes = load(file_linkAttributes) ; %each link results in two states! (two ways)
TableAttributes = readtable(file_linkAttributes2);
TurnAttributes = dlmread(file_turnattributes);

[nbStatesRealStart,nbStatesRealStartAbsorb] = size(incidenceFull); %The states are all possible links
nbRealStates=length(linkAttributes(:,1));
linkNumber = table2cell(TableAttributes(:,1));
linkNumber=cell2mat(linkNumber);
FromNode=table2cell(TableAttributes(:,2));
FromNode=cell2mat(FromNode);
ToNode=table2cell(TableAttributes(:,3));
ToNode=cell2mat(ToNode);
LengthCell=table2cell(TableAttributes(:,5));
LinkLenght=zeros(nbRealStates,1);
for k = 1:nbRealStates
    LinkLength(k,1)=sscanf(LengthCell{k},'%fmi');
end
%conversion from mile to 1/1000 feet 
LinkLength=LinkLength*5.28; 

BikeLane=linkAttributes(:,19);
BikeBlvd=linkAttributes(:,15);
FacType=linkAttributes(:,18);
RMUP=ismember(FacType,0);
BikeFac=BikeLane | BikeBlvd | RMUP;

LargeVolume=linkAttributes(:,9);
HugeVolume1=linkAttributes(:,11);
HugeVolume2=linkAttributes(:,13);
HugeVolume=HugeVolume1 | HugeVolume2;

BridgeNoFac=linkAttributes(:,23); %No bike fac
BridgeNoSepFac=linkAttributes(:,24); %Bike fac NOT separated from traffic, e.g. bike lane
Bridge=linkAttributes(:,25); %All kinds of Bridge
%Bike fac separated from traffic, e.g. PTH-REMU (bike path, rmup) or OTH-SWLK or else
BridgeSepFac=Bridge.*(ones(length(RMUP),1)-BridgeNoFac-BridgeNoSepFac);
BridgeBikeFac= BridgeNoSepFac | BridgeSepFac; 

Upslope4to6p=linkAttributes(:,30);
Upslope6p=linkAttributes(:,32);
Upslope4p=Upslope4to6p | Upslope6p;

%Associate From-Via-To nodes with link pairs:
LinkPairs=zeros(length(TurnAttributes),2);
indices=1:length(linkNumber);
col1=1;
col2=1;
for i=1:length(TurnAttributes)
  z=FromNode==TurnAttributes(i,1);
  y=ToNode==TurnAttributes(i,2);
  t=FromNode==TurnAttributes(i,2);
  x=ToNode==TurnAttributes(i,3);
  if (~isempty(indices(z&y)) && ~isempty(indices(x&t)))
    LinkPairs(col1,1)=indices(z&y);
    LinkPairs(col1,2)=indices(t&x);
    col1=col1+1;
  else
    Todelete(col2,1)=i;
    col2=col2+1;
  end
end
TurnAttributes(Todelete,:)=[]; %some node trios don't correspond to any link pair

LeftTurnMediumTraffic=sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
LeftTurnHighTraffic=sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
NoTurn=sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
LeftTurn=sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
IsLeft=TurnAttributes(:,4)==3;
IsLeftMedTraf=TurnAttributes(:,6) | TurnAttributes(:,10);
IsLeftHighTraf=TurnAttributes(:,7) | TurnAttributes(:,11);
for i=1:length(TurnAttributes)
  LeftTurnMediumTraffic(LinkPairs(i,1),LinkPairs(i,2))=IsLeftMedTraf(i);
  LeftTurnHighTraffic(LinkPairs(i,1),LinkPairs(i,2))=IsLeftHighTraf(i);
  NoTurn(LinkPairs(i,1),LinkPairs(i,2))=1-TurnAttributes(i,9);
  LeftTurn(LinkPairs(i,1),LinkPairs(i,2))=IsLeft(i);
end
  
LinkLength(nbStatesRealStartAbsorb) = 0;
Upslope4p(nbStatesRealStartAbsorb)=0;
LargeVolume(nbStatesRealStartAbsorb)=0;
HugeVolume(nbStatesRealStartAbsorb)=0;
BikeBlvd(nbStatesRealStartAbsorb)=0;
RMUP(nbStatesRealStartAbsorb)=0;
BikeLane(nbStatesRealStartAbsorb)=0;
BridgeBikeFac(nbStatesRealStartAbsorb)=0;
Bridge(nbStatesRealStartAbsorb)=0;

icdLength = incidenceFull * spdiags(LinkLength,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdUpslope4p = incidenceFull * spdiags(Upslope4p,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdLargeVolume = incidenceFull * spdiags(LargeVolume,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdHugeVolume = incidenceFull * spdiags(HugeVolume,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBikeBlvd = incidenceFull * spdiags(BikeBlvd,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdRMUP = incidenceFull * spdiags(RMUP,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBikeLane = incidenceFull * spdiags(BikeLane,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBridgeBikeFac = incidenceFull * spdiags(BridgeBikeFac,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);
icdBridge = incidenceFull * spdiags(Bridge,0,nbStatesRealStartAbsorb,nbStatesRealStartAbsorb);

Atts = objArray(5);
Atts(1).value = icdLength;
Atts(2).value = (icdLargeVolume.*icdLength);
Atts(3).value = (icdHugeVolume.*icdLength);
Atts(4).value = (icdUpslope4p.*icdLength);
Atts(5).value = incidenceFull;
Atts(6).value = (icdBikeBlvd.*icdLength);
Atts(7).value = (icdRMUP.*icdLength);
Atts(8).value = (icdBikeLane.*icdLength);
Atts(9).value = icdBridge;
Atts(10).value = icdBridgeBikeFac;
Atts(11).value = sparse(nbStatesRealStart,nbStatesRealStartAbsorb);
Atts(12).value = LeftTurnMediumTraffic;
Atts(13).value = LeftTurnHighTraffic;
Atts(14).value = NoTurn;
