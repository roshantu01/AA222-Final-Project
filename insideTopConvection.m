function [h] = insideTopConvection(Ts, Ti, r)
%This function calculates the heat transfer coefficient for free convection
%along the inside ceiling
%  Properties are assumed constant and at inside temp because of lack of
%  robust equations evaluating them at a film temperature each time we
%  interate.
L = (pi* r^2)/(2*pi*r);

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
RaL = (g*B*(Ti - Ts)*L^3)/(v*a);
NuL = 0.15*RaL^(1/3);
h = (NuL*k)/L;


end

