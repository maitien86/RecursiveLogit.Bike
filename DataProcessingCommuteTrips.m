%% Data processing
%% Reading the data
dataLinkIDs = xlsread('./Input/LinkIDs.xlsx');
fid = fopen('Observed.csv');
formatSpec = '%s';
N = 38;
dataObs=textscan(fid,repmat('%s',1,38),'HeaderLines',1,'Delimiter',',');
fclose(fid);
LinkID=dataObs{1};
LinkID=str2double(LinkID(1:end));
TripID=dataObs{3};
TripID=str2double(TripID(1:end));
FNode=dataObs{7};
FNode=str2double(FNode(1:end));
TNode=dataObs{8};
TNode=str2double(TNode(1:end));
TripPurpose=dataObs{38};
Ow_rest=dataObs{16};
Ow_rest=str2double(Ow_rest(1:end));
dataObsMat=[LinkID TripID FNode TNode Ow_rest];

%% Link incidence
[nLinks,nbAtts] = size(dataLinkIDs);
temp = dataLinkIDs(:,4);
temp(find(temp == 2)) = 1;
temp = 2 - temp;
% Swap for TF links (one direction)
% for i = 1:nLinks
%     if dataLinkIDs(i,4) == 1
%         dataLinkIDs(i,[2,3]) = dataLinkIDs(i,[3,2]);
%     end
% end

nbiLinks = 2*nLinks;

mappingLinks = zeros(nbiLinks,1); % For mapping
dataFullLinksIDs = zeros(nbiLinks,nbAtts);
dataFullLinksIDs(1:nLinks,:) = dataLinkIDs;
t = nLinks + 1;
for i = 1:nLinks %create the reversed links in any case
    Link = dataLinkIDs(i,:); %link in the original FT direction
    biLink = Link; %We create a TF copy of the link with reversed attributes
    biLink(2) = Link(3);
    biLink(3) = Link(2);
    biLink(7) = Link(8);
    biLink(8) = Link(7);
    biLink(10) = mod(Link(11)+180,360);
    biLink(11) = mod(Link(10)+180,360);       
    if (Link(4)==1) %if there is a FT restriction
        %We switch the order of the TF link (->index i) and FT link (->index t) 
        %The TF link is the one that can be used and the FT link is the
        %restricted one, so we make its utility very high
        dataFullLinksIDs(i,:)=biLink;
        dataFullLinksIDs(t,:)=Link;
        %dataFullLinksIDs(t,6)=10000; %make length absurdly high
    elseif (Link(4) == 2) %if there is a TF restriction
        %we make the utility of the reversed TF link very high
        dataFullLinksIDs(t,:)=biLink;
        %dataFullLinksIDs(t,6)=10000; %make length absurdly high
    else %if Link(4)=0 and there is no restriction
        dataFullLinksIDs(t,:) = biLink;
    end
        mappingLinks(i) = t;
        mappingLinks(t) = i;
        t = t + 1;       
end
%dlmwrite('MappingUTLinksNRL.txt',mappingLinks);

%% connect the network completely
load('unconnectedNodes.mat')
nFakeLinks = length(unconnectedNodes);
for i = 1:nFakeLinks
    %create a certain number of fake links
    %this connect artificially isolated links to the "main" connected
    %component in the network
    Link=dataFullLinksIDs(end,:); %any link can do, we take the last one
    Link(6)=10000; %make length absurdly high
    Link(2)=dataFullLinksIDs(unconnectedNodes(i),3); %the endnode of the unconnected link
    Link(3)=dataFullLinksIDs(10,2); %the startnode of link 10 (which is connected to all destinations!)
    dataFullLinksIDs(nbiLinks+i,:)=Link;       
end
% for i = 1:nFakeLinks
%     %create reverse fake links
%     %this connects the main connected component to the isolated links
%     Link=dataFullLinksIDs(end,:); %any link can do, we take the last one
%     Link(6)=100000; %make length absurdly high
%     Link(2)=dataFullLinksIDs(10,3); %the endnode of link 10
%     Link(3)=dataFullLinksIDs(unconnectedNodes(i),2); %the startnode of the unconnected link
%     dataFullLinksIDs(nbiLinks+nFakeLinks+i,:)=Link;       
% end

%% observations files
FalseObs=[];
tripNum=1;
t = 1;
i = 1;
[sX,sY] = size(dataObs);
sX=length(dataObs{1});
Obs = [];
tripId = [];
while (i<sX)
    j = i;
    while(j~= sX && dataObsMat(j+1,2) == dataObsMat(j,2)) %while it's the same trip ID
        j = j+1; 
    end
    trip = dataObsMat(i:j,:); %lines i to j correspond to a series of links that form a trip
    pathObs = [];
    sizeOfTrip = j-i+1;
    % Take the first sink node
    % Depending on which direction the first link is taken, either the From
    % or the To Node is the sink.
    if ~strcmp(TripPurpose(i),'Exercise')
        if trip(1,4) == trip(2,4) || trip(1,4) == trip(2,3) %if obs goes in FT direction
            sinkNode = trip(1,4);        
            if trip(1,5) == 1; %if FT is not a valid direction for this link (i.e. obs is unvalid)
                pathObs = [pathObs,mappingLinks(trip(1,1))]; %the linkID is the fake link with very high cost
                FalseObs(end+1)=tripNum;
            else
                pathObs = [pathObs,trip(1,1)]; %LinkID is the original linkID in the data
            end
        else %if obs goes in TF direction
            sinkNode = trip(1,3);
            if trip(1,5) == 1 %if FT is not a valid direction
                pathObs = [pathObs,trip(1,1)]; %LinkID corresponding to TF is the original linkID in the data            
            else %if FT is a valid direction
                pathObs = [pathObs,mappingLinks(trip(1,1))]; %LinkID corresponding to TF is the matching reversed link
                if trip(1,5)==2
                    FalseObs(end+1)=tripNum;
                end
            end
        end
        for k = 2:sizeOfTrip %do that again for all the remaining links in the trip
           if (trip(k,3) == sinkNode)
                sinkNode = trip(k,4);
                if trip(k,5) ==1;
                    pathObs = [pathObs,mappingLinks(trip(k,1))];
                    FalseObs(end+1)=tripNum;
                else
                    pathObs = [pathObs,trip(k,1)];
                end
           else
                sinkNode = trip(k,3);
                if trip(k,5) == 1;
                    pathObs = [pathObs,trip(k,1)];                
                else
                    pathObs = [pathObs,mappingLinks(trip(k,1))]; 
                    if trip(k,5) == 2
                        FalseObs(end+1) = tripNum;
                    end
                end
           end
        end
        pathObs  = [pathObs(end), pathObs];%a trick useful later
        Obs(end+1,1:size(pathObs,2)) = pathObs;
        tripId = [tripId, dataObsMat(i,2)];
        tripNum=tripNum+1;
    end
    i = j+1;       
end
Obs = sparse(Obs);
% To remove false obs - but then not enough data for NRL
% uniqueFalseObs=unique(FalseObs); 
% Obs(uniqueFalseObs,:)=[];

%% Creating incidence matrix ###
temp = find(Obs(:,1)); %this corresponds to the destination because we added pathObs(end) to the beginning of pathObs!
Dest = unique(Obs(temp,1));
nDummyLinks = size(Dest,1);
MappingDestDummy = zeros(1,nbiLinks);
A = dataFullLinksIDs(:,3); %To Nodes (end of link)
B = dataFullLinksIDs(:,2); %From Nodes (start of link)
Incidence = sparse(nbiLinks+nFakeLinks,nbiLinks+nFakeLinks + nDummyLinks);
for i = 1:nbiLinks+nFakeLinks
    U = find(B == A(i));%for each link i, look at all the links U that connect to the End Node.
    Incidence(i,U) = 1;% this means we can go from link i to link U
end
for i = 1: nDummyLinks
    Incidence(Dest(i),nbiLinks+nFakeLinks + i) = 1;
    MappingDestDummy(Dest(i)) = nbiLinks+nFakeLinks + i;
end

%% Remapping Obs 
temp = Obs(:,1);
temp = MappingDestDummy(temp)';
Obs(:,1) = temp; %adding the fake absorbing link after each observed path
for i =1: size(Obs,1)
    d = size(find(Obs(i,:)),2);
    Obs(i,d+1) = temp(i);
end

%% Save data to files

[i,j,val] = find(Obs);
data_dump = [i,j,val];
save('TripsObservationsNRLNonExercise.txt','data_dump','-ascii');

[i,j,val] = find(Incidence);
data_dump = [i,j,val];
save('IncidenceNRLNonExercise.txt','data_dump','-ascii');

dlmwrite('fullLinksAttsNRLNonExercise.txt',dataFullLinksIDs);

[i,j,val] = find(MappingDestDummy);
data_dump = [i,j,val];
save('mappingDestDummyNRLNonExercise.txt','data_dump','-ascii');

