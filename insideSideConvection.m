function [h] = insideSideConvection(Ts, Ti, h)
%This function calculates the heat transfer coefficient for free convection
%along the inside walls
%  Properties are assumed constant and at inside temp because of lack of
%  robust equations evaluating them at a film temperature each time we
%  interate.

L = h;
Pr = 0.713;
dens = 0.818;
u = 18.5*10^-6;
v = u/dens;
Cp = 1010;
k = 0.0261;
a = k/(dens*Cp);
g = 3.7;
B = 0.00336;

%RaL number eq and Nu correlation from Appendix A
RaL = ((g*B*(Ti - Ts)*L^3)/(v^2)) * Pr;
NuL = ( 0.825 + ( (0.387*RaL^(1/6))/((1+(0.492/Pr)^(9/16))^(8/27))))^2;

h = (NuL*k)/L;



end

