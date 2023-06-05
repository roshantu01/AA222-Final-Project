function [hr] = outsideRadiation(Ts,Tsurr, E)
%This function calculates the effective heat tranfer coefficient due to
%radiation, in accordance with the eq from Appendix A
% E = 0.03;
sigma = 5.67*10^-8;

hr = E*sigma*(Ts + Tsurr) * (Ts^2 + Tsurr^2);


end

