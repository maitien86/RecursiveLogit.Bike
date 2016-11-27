%% Data processing
%% Reading the data
dataLinkIDs = xlsread('./Input/BikeFacility.xlsx');

%% Link incidence
[nLinks,nbAtts] = size(dataLinkIDs);
nbiLinks = 2*nLinks;
mappingLinks = zeros(nbiLinks,1); % For mapping
%write in original format
dataFullLinksBis = zeros(nLinks,nbAtts);
dataFullLinksIDs = dataLinkIDs;
dataFullLinksIDs(:,10:11)=dataLinkIDs(:,4:5); %fbearing,tbearing
dataFullLinksIDs(:,7:8)=dataLinkIDs(:,18:19); %dzp dzn
dataFullLinksIDs(:,6)=dataLinkIDs(:,12); %length
dataFullLinksIDs(:,5)=dataLinkIDs(:,13); %ow_rest
dataFullLinksIDs(:,26)=dataLinkIDs(:,21); %est vol
dataFullLinksIDs(:,33)=dataLinkIDs(:,28); %bridge
dataFullLinksIDs=vertcat(dataFullLinksIDs,dataFullLinksBis);
t = nLinks + 1;
for i = 1:nLinks %create the reversed links in any case
    Link = dataFullLinksIDs(i,:); %link in the original FT direction
    biLink = Link; %We create a TF copy of the link with reversed attributes
    biLink(2) = Link(3);
    biLink(3) = Link(2);
    biLink(7) = Link(8);
    biLink(8) = Link(7);
    biLink(10) = mod(Link(11)+180,360);
    biLink(11) = mod(Link(10)+180,360);       
    if (strcmp(Link(5),'TF')) %if there is a FT restriction
        %We switch the order of the TF link (->index i) and FT link (->index t) 
        %The TF link is the one that can be used and the FT link is the
        %restricted one, so we make its utility very high
        dataFullLinksIDs(i,:)=biLink;
        dataFullLinksIDs(t,:)=Link;
        %dataFullLinksIDs(t,6)=10000; %make length absurdly high
    elseif (strcmp(Link(5),'FT')) %if there is a TF restriction
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

%% Find at which nodes to add dummy links
filename='./Input/TAZ_Node_Link.xlsx';
Table=readtable(filename);
CentroidToNode=table2array(Table(:,1:2));
CentroidToNode(end+1,:)=[24,0]; %row 24 missing
CentroidToNode2=sortrows(CentroidToNode);
NbOD=length(CentroidToNode(:,1));

% Match each OD node to links
DestLinks=cell(1,NbOD);
OrigLinks=cell(1,NbOD);
for k=1:NbOD
    DestLinks{k}=find(dataFullLinksIDs(:,3)==CentroidToNode2(k,2));
    OrigLinks{k}=find(dataFullLinksIDs(:,2)==CentroidToNode2(k,2));
end

%% Creating incidence matrix ###
nStartLinks = double(NbOD);
nAbsorbingLinks = double(NbOD);
A = dataFullLinksIDs(:,3); %To Nodes (end of link)
B = dataFullLinksIDs(:,2); %From Nodes (start of link)
Incidence = sparse(nbiLinks+nStartLinks,nbiLinks+nStartLinks+nAbsorbingLinks);
for i = 1:nbiLinks
    U = find(B == A(i));%for each link i, look at all the links U that connect to the End Node.
    Incidence(i,U) = 1;% this means we can go from link i to link U
end
for i = 1: nAbsorbingLinks
    Incidence(DestLinks{i},nbiLinks + nStartLinks + i) = 1;
end
for i = 1: nStartLinks
    Incidence(nbiLinks + i,OrigLinks{i}) = 1;
end

%% Read OD demand matrix
filename='./Input/bucket_rounded_bike_new.txt';
fid=fopen(filename);
[Nout]=textscan(fid,'%d%d:%d%d:%d%d:%d%d:%d%d:%d','CollectOutput',1,'EmptyValue',0);
fclose(fid);
out=[Nout{:}];
NbOD=out(end,1);
ODdemand=zeros(NbOD);
%Shape into demand matrix indexed with with centroid ID
for i=1:length(out(:,1))
    for j=1:5
        if out(i,2*j) ~= 0
            ODdemand(out(i,1),out(i,2*j))=out(1,2*j+1);
        end
    end
end

%% Write OD demand in terms of link IDs
[row, col]=find(ODdemand>0);
ODmat=sparse(zeros(size(Incidence)));
for k=1:length(row)
    i=row(k);%origin Centroid ID
    j=col(k);%destination Centroid ID
    ODmat(nbiLinks+i,nbiLinks+nStartLinks+j)=ODdemand(i,j);
end

%% Write data to file
[i,j,val] = find(Incidence);
data_dump = [i,j,val];
save('IncidenceForLinkFLows.txt','data_dump','-ascii');

[i,j,val] = find(ODmat);
data_dump = [i,j,val];
save('ODmat.txt','data_dump','-ascii');

dlmwrite('fullLinksAttsForLinkFlows.txt',dataFullLinksIDs);

[i,j,val] = find(mappingLinks);
data_dump = [i,j,val];
save('mappingLinksForLinkFlows.txt','data_dump','-ascii');

