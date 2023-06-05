function [mass] = mass_calc(x)
    r = 10;
    h = 8;
    %Lengths
    length_MLI = x(1);
    length_rego = x(2);
    length_aero = x(3);
    
    length_tot = sum(x);
    
    %Densities
    p_MLI = 60;
    p_aero = 130;
    p_rego = 1800;

    %side/tube Volumes
    vs_aero = tubevol(r, r+length_aero,h);
    vs_rego = tubevol(r+length_aero, r+length_aero+length_rego,h);
    vs_MLI = tubevol(r+length_aero+length_rego, r + length_tot,h);

    %top volumes
    vt_aero = pi*((r+length_tot)^2)*length_aero;
    vt_rego = pi*((r+length_tot)^2)*length_rego;
    vt_MLI = pi*((r+length_tot)^2)*length_MLI;

    mass = (vs_aero + vt_aero)*p_aero + (vs_rego + vt_rego)*p_rego + (vs_MLI + vt_MLI)*p_MLI;



end

function v = tubevol(ri, ro, h)
    v = (pi*ro^2 - pi*ri^2)*h;
end