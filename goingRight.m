function [ RightTurn ] = goingRight( angle1,angle2 )
% %Going Right -- FALSE, right and left need to be switched
% %   Determines whether the bike makes a right turn given the 2 angles.
% if (angle1-180 >=0) %if angle1 is greater (or equal) than 180
%     RightTurn=(angle1-angle2>30) && (angle1-180-angle2<0); %angle 2 is within (angle1-180,angle1-30)
% elseif (angle1-30>0) %if angle1 is smaller than 180 but also greater than 30
%     RightTurn=(angle2-angle1-180>0) || (angle1-angle2>30); %angle 2 is either within (angle1+180,360] or within [0,angle1-30)
% else %if angle1 is smaller than 30
%     RightTurn=(angle2-angle1-180>0) && (angle2-angle1-330<0) ;%angle2 is within (angle1+180,angle1+330)
% end

%Going Left 
%   Determines if the bike makes a left turn given the angle of the
%   outgoing link (angle1) and the incoming link (angle2)
%   Return 1 if there is a left turn, 0 else.
if (angle1 <180) %if angle1 is smaller than 180
    RightTurn=(angle2-angle1>60) && (angle1+180-angle2>0); %angle 2 is within (angle1+60,angle1+180)
elseif (angle1 < 300) %if angle1 is greater (or equal) than 180 but smaller than 60
    RightTurn= (angle2-angle1-60>0) || (angle1-180-angle2>0);%angle 2 is either within (angle1+60,360] or within [0,angle1-180)
else %if angle1 is between 300 and 360 (included)
    RightTurn=(angle2-angle1+300>0)&&(angle1-angle2-180>0);%angle2 is within (angle1-300,angle1-180)
end

end

