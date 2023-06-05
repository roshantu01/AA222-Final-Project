function [q] = heat_loss(x)
    length_MLI = x(1);
    length_aero = x(3);
    length_rego = x(2);
    length_tot = sum(x);
    
    %% Defining constants and some initial guesses
    Tinf_outside = 148.15;
    Tsidesurf_outside = 200; %initial guess
    Ttopsurf_outside = 200; %initial guess
    
    pressure = 600;
    velocity = 30;
    
    %Surface Areas used for computation
    r_inner = 10;
    r_outer = r_inner+length_tot;
    h = 8;

    A_side_outer = 2*pi*r_outer*h;
    A_side_inner = 2*pi*r_inner*h;

    A_top_outer = pi*r_inner^2;
    A_top_inner = pi*r_inner^2;
    %Defining heat transfer params based on x
    E = 0.03/(length_MLI*1000);
    k_aero = 0.014;
    k_rego = 0.4;

    
    %% Outside sides Convection
    [hos, Pr] = outerSideConvection(Tinf_outside, Tsidesurf_outside, pressure, velocity, r_outer);
    %Defining thermal resistance
    Rcos = 1/(hos*A_side_outer);
    %% Outside Top Convection
    hot = outerTopConvection(Tinf_outside, Ttopsurf_outside, pressure, velocity, r_outer);
    %Defining thermal resistance
    Rcot = 1/(hot*A_top_outer);
    %% Inside sides Convection
    Tsidesurf_inside = 250; %initial guess
    Tinf_inside = 298.15;
    his = insideSideConvection(Tsidesurf_inside, Tinf_inside, h);
    %Defining thermal resistance
    Rcis = 1/(his*A_side_inner);
    
    %% Inside Top Convection
    Ttopsurf_inside = 250;
    hit = insideTopConvection(Ttopsurf_inside, Tinf_inside, r_inner);
    %Defining thermal resistance
    Rcit = 1/(hit*A_top_inner);
    %% Outside Radiation
    hros = outsideRadiation(Tsidesurf_outside, Tinf_outside, E);
    hrot = outsideRadiation(Ttopsurf_outside, Tinf_outside, E);
    %Defining thermal resistances
    Rros = 1/(hros*A_side_outer);
    Rrot = 1/(hrot*A_top_inner);
    %% Conduction
%     k = 0.014;
    %Defining thermal resistance for aerogel
    Rconds_aero = log((r_inner+length_aero)/r_inner)/(2*pi*k_aero*h);
    Rcondt_aero = length_aero/(k_aero*0.5*(A_top_outer+A_top_inner));

    Rconds_rego = log((r_inner+length_aero+length_rego)/(r_inner+length_aero))/(2*pi*k_rego*h);
    Rcondt_rego = length_rego/(k_rego*0.5*(A_top_outer+A_top_inner));

    Rconds = Rconds_aero+Rconds_rego;
    Rcondt = Rcondt_aero+Rcondt_rego;


    %% First iteration using intial guesses
    R1 = ((1/Rros) + (1/Rcos))^-1;
    Rtotside = Rcis + Rconds + R1;
    Qside = (Tinf_inside - Tinf_outside)/Rtotside;
    R2 = ((1/Rrot) + (1/Rcot))^-1;
    Rtottop = Rcit + Rcondt + R2;
    Qtop = (Tinf_inside - Tinf_outside)/Rtottop;
    Qtot = Qside+Qtop;
    
    Qtotinitial = Qtot;
    
    Qset = [Qtot];
    Tsidesurf_insideset = [Tsidesurf_inside];
    Ttopsurf_insideset = [Ttopsurf_inside];
    Tsidesurf_outsideset = [Tsidesurf_outside];
    Ttopsurf_outsideset = [Ttopsurf_outside];
    iterations = [0];
    
    %% Iterating solution, updating itermediate temps
    for i = 1:10
        %inside side surface temp using inside convection resistance
        Tsidesurf_inside = Tinf_inside - (Qside*Rcis);
        his = insideSideConvection(Tsidesurf_inside, Tinf_inside, h);
        Rcis = 1/(his*A_side_inner);
        %outside side surface temp using outside convection and radiation resistance
        Tsidesurf_outside = Qside*R1 + Tinf_outside;
        [hos, Pr] = outerSideConvection(Tinf_outside, Tsidesurf_outside, pressure, velocity, r_outer);
        Rcos = 1/(hos*A_side_outer);
        hros = outsideRadiation(Tsidesurf_outside, Tinf_outside, E);
        Rros = 1/(hros*A_side_outer);
    
        R1 = ((1/Rros) + (1/Rcos))^-1;
        Rtotside = Rcis + Rconds + R1;
        Qside = (Tinf_inside - Tinf_outside)/Rtotside;
        
        
        
        %inside top surface temp using inside convection resistance
        Ttopsurf_inside = Tinf_inside - (Qtop*Rcit);
        hit = insideTopConvection(Ttopsurf_inside, Tinf_inside, r_inner);
        Rcit = 1/(hit*A_top_inner);
        %outside top surface temp updated using outside convection and
        %radiation
        Ttopsurf_outside = Qtop*R2 + Tinf_outside;
        hot = outerTopConvection(Tinf_outside, Ttopsurf_outside, pressure, velocity, r_outer);
        Rcot = 1/(hot*A_top_outer);
        hrot = outsideRadiation(Ttopsurf_outside, Tinf_outside, E);
        Rrot = 1/(hrot*A_top_outer);
        
        R2 = ((1/Rrot) + (1/Rcot))^-1;
        Rtottop = Rcit + Rcondt + R2;
        Qtop = (Tinf_inside - Tinf_outside)/Rtottop;
        Qtot = Qside+Qtop;
        
        Qset = [Qset Qtot];
        Tsidesurf_insideset = [Tsidesurf_insideset Tsidesurf_inside];
	    Ttopsurf_insideset = [Ttopsurf_insideset Ttopsurf_inside];
        Tsidesurf_outsideset = [Tsidesurf_outsideset Tsidesurf_outside];
        Ttopsurf_outsideset = [Ttopsurf_outsideset Ttopsurf_outside];
        iterations = [iterations i];
        
        
    end
    
    %% Plotting Values vs Iterations to show convergence;
    q = Qtot;

% figure;
% hold on
% plot(iterations, Qset);
% plot(iterations, Qset, 'ro');
% xlabel("Iterations");
% ylabel("Total Heat Loss [W]");
% title("Total Heat Loss");
% 
% 
% figure;
% hold on;
% plot(iterations, Tsidesurf_insideset);
% plot(iterations, Ttopsurf_insideset);
% plot(iterations, Tsidesurf_outsideset);
% plot(iterations, Ttopsurf_outsideset);
% legend("Inside Wall", "Inside Ceiling", "Outside Wall", "Outside Ceiling");
% xlabel("Iterations");
% ylabel("Temperature [K]");
% title("Surface Temperatures");

end











