function [h] = outerTopConvection(To, Ts, p, Vel,r)
%OuterTopConvection calculates and returns the convection coefficient of
%the outside top of the habitat
%   Atmospheric properties are calculated using equations in Appendix A
    D = (pi* r^2)/(2*pi*r);

    Tf = (To+Ts)/2;
    
    dens = p/(188.918 * Tf);
    u = (1.37*10^-5)* ((Tf/273)^1.5) * ( (273+222)/(Tf+222) );
    v = u/dens;
    k = -7.215*10^-3 + ((8.015*10^-5)*(Tf)) + ((5.477*10^-9)*(Tf^2)) + ((-1.053*10^-11)*(Tf^3));
    Cp2 = (4.34247158*10^2) + (1.69739101*Tf) + ((-1.34580010*10^-3)*(Tf^2)) + ((4.64595961*10^-7)*(Tf^3)) + ((-2.71480543*10^-11)*(Tf^4));

    %Pr and Re calculated using these values
    Pr = (Cp2*u)/k;
    Re = (Vel*D)/v;
    
    %Nusselt number correlation from Appendix A
    NuL = (0.037*Re^(4/5))*Pr^(1/3);
    h = (NuL*k)/D;
    
end

