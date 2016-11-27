%   Compute Link Size attribute from data
%   
%%
function ok = getGlobalLinkFlow()
    global incidenceFull; 
    global LSatt;
    global GLS;
    global nbobs;
    global EstimatedTime;
    
    importfile('../LSatt.mat');
    GLS = sparse(zeros(size(incidenceFull)));
    for n = 1:nbobs
        n
        GLS = sparse(GLS + LSatt(n).value);
    end    
    GLS = GLS/1000;
    ok = true;
end