function [ noTurn ] = goingStraight( angle1, angle2 )
%goingStraight 
%   Determines whether the bike goes straight without turning (±5 degrees)
%   given angles of outgoing and incoming links (angle1 and angle2)
if (angle1+5 <360 && angle1-5>0) %if angle1 is within (5,355)
    noTurn=(angle2-angle1<5) && (angle2-angle1>-5); %angle 2 is within (angle1-5,angle1+5)
elseif (angle1+5 >= 360) %if angle1 is greater than 355
    noTurn=(angle1-355-angle2>0)||(angle2>angle1-5); %angle2 is within 5deg of angle 1
else %if angle1 is smaller than 5
    noTurn=(angle2>angle1+355) || (angle2<angle1+5); %angle2 is within 5deg of angle 1
end

end

