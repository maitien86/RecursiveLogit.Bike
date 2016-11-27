function [] = CreatSubNetWork(nRemovedLink)
    global incidenceFull;
    global Atts;  
    nbLinks = size(incidenceFull,1);
    I = find(incidenceFull);
    nbnonzero = size(I,1);
    %newIcd = incidenceFull;
    newObs = [];
    %% Removing links
    rM = randi([1,nbnonzero],1,nRemovedLink);
    t = 0;
    for i=1: nRemovedLink
        [k,a] = ind2sub(size(incidenceFull), I(rM(i)));
        if a <= 7288 && size(find(incidenceFull(k,:)),2) >1  
            incidenceFull(k,a) = 0;
            t = t+1;
            Atts = getAtt();
            if (getLL_checkConnectedNetwork() == false)
                incidenceFull(k,a) = 1;
                t = t-1;
            end
        end
        t
    end
    
    %incidenceFull = newIcd;
    
    %getLL_checkConnectedNetwork();
    %% Resample sample
%     nbobs = size(Observations,1)
%     for i=1:nbobs
%         path = Observations(i,:);
%         ok = true;
%         for j=2: size(find(path),2)-1
%             if newIcd(path(j), path(j+1)) == 0
%                 ok = false;
%             end
%         end
%         if ok == true
%             newObs = [newObs;path];
%         end
%     end
end


