%   Get attribute 
function Att = getAtt()

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;
    global GLS;
    global isLinkSizeInclusive;
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    mu = 1;
    Incidence = 1/mu * incidenceFull;
    Att(1) = Matrix2D(Incidence .* EstimatedTime);
  %  Att(1).Value = Att(1).Value(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Att(2) = Matrix2D(Incidence .* TurnAngles);
 %   Att(2).Value = Att(2).Value(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Att(3) = Matrix2D(Incidence .* LeftTurn);
 %   Att(3).Value = Att(3).Value(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Att(4) = Matrix2D(Incidence .* Uturn);
 %   Att(4).Value = Att(4).Value(1:lastIndexNetworkState,1:lastIndexNetworkState);
    if isLinkSizeInclusive == true
        Att(5) = Matrix2D(Incidence .* GLS);
    end
end
