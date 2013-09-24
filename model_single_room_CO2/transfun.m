function [ Crdot ] = transfun( t,Cr )

Ci = 1; %input CO2 conc.
V = 100; % volume of room
Ca = 2; %ambient CO2 of outside office
 

if ( t >= 10 & t <= 25 )
    Fo = 0;
else
    Fo = 10; % output flow in m^3/s
end
if ( t > 15 & t <= 40 )
    Fi = 10;
else
    Fi = 0;
end
Fdoor = Fi - Fo;


if ( Fdoor > 0 )
   doorterm = -Fdoor*Cr; 
else
   doorterm = -Fdoor*Ca;
end


Crdot = 1/V * (-Fo*Cr + doorterm + Ci*Fi + 1);
end

