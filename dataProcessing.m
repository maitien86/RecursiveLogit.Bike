%% Data processing
% Creating Data Files
% Original Network with restrictions in given directions
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

nbiLinks = sum(temp); %nbiLinks is the total number of links including return links but not fake absorbing links

mappingLinks = zeros(nbiLinks,1); % For mapping
dataFullLinksIDs = zeros(nbiLinks,nbAtts);
dataFullLinksIDs(1:nLinks,:) = dataLinkIDs;
t = nLinks + 1;
for i = 1:nLinks %create the reverse links if they exist
    Link = dataLinkIDs(i,:);
    if (Link(4) == 0)
        biLink = Link;
        biLink(2) = Link(3);
        biLink(3) = Link(2);
        biLink(7) = Link(8);
        biLink(8) = Link(7);
        biLink(10) = mod(Link(11)+180,360);
        biLink(11) = mod(Link(10)+180,360);
        dataFullLinksIDs(t,:) = biLink;
        mappingLinks(i) = t;
        mappingLinks(t) = i;
        t = t + 1;
    end     
end
for i = 1:nLinks
    if dataLinkIDs(i,4)==1 %if there is a FT restriction, reverse attributes
        tempo=dataFullLinksIDs(i,2);
        dataFullLinksIDs(i,2)=dataFullLinksIDs(i,3);
        dataFullLinksIDs(i,3)=tempo;
        tempo=dataFullLinksIDs(i,7);
        dataFullLinksIDs(i,7)=dataFullLinksIDs(i,8);
        dataFullLinksIDs(i,8)=tempo;
        tempo=dataFullLinksIDs(i,10);
        dataFullLinksIDs(i,10)=mod(dataFullLinksIDs(i,11)+180,360);
        dataFullLinksIDs(i,11)=mod(tempo+180,360);
    end
end
%dlmwrite('MappingUTLinks.txt',mappingLinks);

%% observations files
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
    if trip(1,4) == trip(2,4) || trip(1,4) == trip(2,3) %if obs goes FT
        sinkNode = trip(1,4); 
        pathObs = [pathObs,trip(1,1)]; %LinkID for FT is the original linkID in the data
        %NOTE: we did not test whether FT is a valid direction for this
        %link (obs might be unvalid). If FT does not exist trip(1,1) will
        %correspond to the ID of the TF link and the resulted path will be
        %unconnected -> need at the end to check whether obs are valid
    else %if obs goes TF
        sinkNode = trip(1,3);
        if trip(1,5) == 0 %if there is no one way restr.
            pathObs = [pathObs,mappingLinks(trip(1,1))]; %LinkID for TF is the reversed "matching" linkID in the copied data
        else %if there is a one way restr.
            pathObs = [pathObs,trip(1,1)];%LinkID is the original linkID which corresponds to TF in this case
            dataLinkIDs(trip(1,1),[2,3]) = dataLinkIDs(trip(1,1),[3,2]);
            dataLinkIDs(trip(1,1),[7,8]) = dataLinkIDs(trip(1,1),[8,7]);           
        end
        %NOTE: we did not check whether TF is a valid direction for this
        %link (obs might be unvalid). If TF does not exist trip(1,1) still
        %corresponds to the ID of the FT link and the resulted path will be
        %unconnected -> need at the end to check whether obs are valid
    end
    for k = 2:sizeOfTrip %do that again for all the remaining links in the trip
       if (trip(k,3) == sinkNode)
            sinkNode = trip(k,4);
            pathObs = [pathObs,trip(k,1)];
       else
            sinkNode = trip(k,3);
            if trip(k,5) == 0
                pathObs = [pathObs,mappingLinks(trip(k,1))];
            else
                pathObs = [pathObs,trip(k,1)];
                dataLinkIDs(trip(k,1),[2,3]) = dataLinkIDs(trip(k,1),[3,2]);
                dataLinkIDs(trip(k,1),[7,8]) = dataLinkIDs(trip(k,1),[8,7]);           
            end
       end
    end
    pathObs  = [pathObs(end), pathObs];%a trick useful later
    Obs(end+1,1:size(pathObs,2)) = pathObs;
    tripId = [tripId, dataObsMat(i,2)];
    i = j+1;    
end
Obs = sparse(Obs);

%% Creating incidence matrix ###
temp = find(Obs(:,1)); %this corresponds to the destination because we added pathObs(end) to the beginning of pathObs!
Dest = unique(Obs(temp,1));
nDummyLinks = size(Dest,1);
MappingDestDummy = zeros(1,nbiLinks);
A = dataFullLinksIDs(:,3); %To Nodes (end of link)
B = dataFullLinksIDs(:,2); %From Nodes (start of link)
Incidence = sparse(nbiLinks,nbiLinks + nDummyLinks);
for i = 1:nbiLinks
    U = find(B == A(i));%for each link i, look at all the links U that connect to the End Node.
    Incidence(i,U) = 1;% this means we can go from link i to link U
end
for i = 1: nDummyLinks
    Incidence(Dest(i),nbiLinks + i) = 1;
    MappingDestDummy(Dest(i)) = nbiLinks + i;
end

%% Remappping Obs 
temp = Obs(:,1);
temp = MappingDestDummy(temp)';
Obs(:,1) = temp; %adding the fake absorbing link after each observed path
for i =1: size(Obs,1)
    d = size(find(Obs(i,:)),2);
    Obs(i,d+1) = temp(i);
end

%% Save data to files

% [i,j,val] = find(Obs);
% data_dump = [i,j,val];
% save('TripsObservations.txt','data_dump','-ascii');
% 
% 
% [i,j,val] = find(Incidence);
% data_dump = [i,j,val];
% save('Incidence.txt','data_dump','-ascii');
% 
% 
% dlmwrite('fullLinksAtts.txt',dataFullLinksIDs);
% 
% [i,j,val] = find(MappingDestDummy);
% data_dump = [i,j,val];
% save('mappingDestDummy.txt','data_dump','-ascii');


