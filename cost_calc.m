function [cost] = cost_calc(x)
    r = 10;
    h = 8;
    %Lengths
    length_MLI = x(1);
    length_rego = x(2);
    length_aero = x(3);
    
    length_tot = sum(x);
    
    %Densities
    c_MLI = 1203000;
    c_aero = 2690416.886;
    c_rego = 1398.733333;

    %side/tube Volumes
    vs_aero = tubevol(r, r+length_aero,h);
    vs_rego = tubevol(r+length_aero, r+length_aero+length_rego,h);
    vs_MLI = tubevol(r+length_aero+length_rego, r + length_tot,h);

    %top volumes
    vt_aero = pi*((r+length_tot)^2)*length_aero;
    vt_rego = pi*((r+length_tot)^2)*length_rego;
    vt_MLI = pi*((r+length_tot)^2)*length_MLI;

    cost = (vs_aero + vt_aero)*c_aero + (vs_rego + vt_rego)*c_rego + (vs_MLI + vt_MLI)*c_MLI;



end

function v = tubevol(ri, ro, h)
    v = (pi*ro^2 - pi*ri^2)*h;
end