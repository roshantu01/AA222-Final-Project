%% Defining constants and some initial guesses
Tinf_outside = 148.15;
Tsidesurf_outside = 200; %initial guess
Ttopsurf_outside = 200; %initial guess

pressure = 600;
velocity = 30;

%Surface Areas used for computation
A_side = 2*pi*10*8;
A_top = pi*10^2;

%% Outside sides Convection
[hos, Pr] = outerSideConvection(Tinf_outside, Tsidesurf_outside, pressure, velocity);
%Defining thermal resistance
Rcos = 1/(hos*A_side);
%% Outside Top Convection
hot = outerTopConvection(Tinf_outside, Ttopsurf_outside, pressure, velocity);
%Defining thermal resistance
Rcot = 1/(hot*A_top);
%% Inside sides Convection
Tsidesurf_inside = 250; %initial guess
Tinf_inside = 298.15;
his = insideSideConvection(Tsidesurf_inside, Tinf_inside);
%Defining thermal resistance
Rcis = 1/(his*A_side);

%% Inside Top Convection
Ttopsurf_inside = 250;
hit = insideTopConvection(Ttopsurf_inside, Tinf_inside);
%Defining thermal resistance
Rcit = 1/(hit*A_top);
%% Outside Radiation
hros = outsideRadiation(Tsidesurf_outside, Tinf_outside);
hrot = outsideRadiation(Ttopsurf_outside, Tinf_outside);
%Defining thermal resistances
Rros = 1/(hros*A_side);
Rrot = 1/(hrot*A_top);
%% Conduction
k = 0.014;
%Defining thermal resistance
Rconds = log(10/9.5)/(2*pi*k*8);
Rcondt = 0.5/(k*A_top);
%% First iteration using intial guesses
R1 = ((1/Rros) + (1/Rcos))^-1;
Rtotside = Rcis + Rconds + R1;
Qside = (Tinf_inside - Tinf_outside)/Rtotside;
R2 = ((1/Rrot) + (1/Rcot))^-1;
Rtottop = Rcit + Rcondt + R2;
Qtop = (Tinf_inside - Tinf_outside)/Rtottop;
Qtot = Qside+Qtop;

Qtotinitial = Qtot

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
    his = insideSideConvection(Tsidesurf_inside, Tinf_inside);
    Rcis = 1/(his*A_side);
    %outside side surface temp using outside convection and radiation resistance
    Tsidesurf_outside = Qside*R1 + Tinf_outside;
    [hos, Pr] = outerSideConvection(Tinf_outside, Tsidesurf_outside, pressure, velocity);
    Rcos = 1/(hos*A_side);
    hros = outsideRadiation(Tsidesurf_outside, Tinf_outside);
    Rros = 1/(hros*A_side);

    R1 = ((1/Rros) + (1/Rcos))^-1;
    Rtotside = Rcis + Rconds + R1;
    Qside = (Tinf_inside - Tinf_outside)/Rtotside;
    
    
    
    %inside top surface temp using inside convection resistance
    Ttopsurf_inside = Tinf_inside - (Qtop*Rcit);
    hit = insideTopConvection(Ttopsurf_inside, Tinf_inside);
    Rcit = 1/(hit*A_top);
    %outside top surface temp updated using outside convection and
    %radiation
    Ttopsurf_outside = Qtop*R2 + Tinf_outside;
    hot = outerTopConvection(Tinf_outside, Ttopsurf_outside, pressure, velocity);
    Rcot = 1/(hot*A_top);
    hrot = outsideRadiation(Ttopsurf_outside, Tinf_outside);
    Rrot = 1/(hrot*A_top);
    
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
Qtot

figure;
hold on
plot(iterations, Qset);
plot(iterations, Qset, 'ro');
xlabel("Iterations");
ylabel("Total Heat Loss [W]");
title("Total Heat Loss");


figure;
hold on;
plot(iterations, Tsidesurf_insideset);
plot(iterations, Ttopsurf_insideset);
plot(iterations, Tsidesurf_outsideset);
plot(iterations, Ttopsurf_outsideset);
legend("Inside Wall", "Inside Ceiling", "Outside Wall", "Outside Ceiling");
xlabel("Iterations");
ylabel("Temperature [K]");
title("Surface Temperatures");













