function [h, Pr] = outerSideConvection(To, Ts, p, Vel, r)

%OuterSideConvection calculates and returns the convection coefficient of
%the outside sides of the habitat
%   Atmospheric properties are calculated using equations in Appendix A
    D = 2*r;
    Tf = (To+Ts)/2;
    
    dens = p/(188.918 * Tf);
    u = (1.37*10^-5)* ((Tf/273)^1.5) * ( (273+222)/(Tf+222) );
    v = u/dens;
    k = -7.215*10^-3 + ((8.015*10^-5)*(Tf)) + ((5.477*10^-9)*(Tf^2)) + ((-1.053*10^-11)*(Tf^3));
    Cp2 = (4.34247158*10^2) + (1.69739101*Tf) + ((-1.34580010*10^-3)*(Tf^2)) + ((4.64595961*10^-7)*(Tf^3)) + ((-2.71480543*10^-11)*(Tf^4));

    %Pr and Re calculated using these values
    Pr = (Cp2*u)/k;
    Re = (Vel*D)/v;

    NuD = 0.3 + ( (0.62*(Re^0.5)*(Pr^(1/3))) / (1+((0.4/Pr)^(2/3)))^0.25 ) * ((1+((Re/282000)^(5/8)))^(4/5)) ;
    
    %Nusselt number correlation from Appendix A
    hc = (NuD*k)/D;
    h = hc;
end

