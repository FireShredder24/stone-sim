% Blowdown rocket system simulator
% (c) John Nguyen 2025
% MIT License

% This program assumes that, within reasonable O/F ratio bounds,
% the rocket chamber pressure is proportional to the mass flow rate
% and the thrust is proportional to the chamber pressure

% If "isothermal" is True, the expansion of the ullage gas in the propellant tanks
% tanks is assumed to be isothermal (constant temperature).
% This assumption is only valid when the expansion cooling of the gas
% is insignificant compared to the heat transfer from the environment
% as may be the case in a poorly insulated static fire test stand
% which exhausts its propellant slowly over a long burn.

% If "isothermal" is False, the gas expansion is assumed to be adiabatic,
% (no heat transfer).  This will result in reduced pressures compared
% to the isothermal case, and is nearly a "worst-case" scenario.  However
% there is still the case of "negative" heat transfer to the gas, which
% can easily occur if the ullage gas was initially at a higher temperature
% than the cryogenic propellant (say, liquid hydrogen) it pressurized.

isothermal = false;

% Common specific heat ratios:
% Helium - 1.667
% Hydrogen - 1.41
% Nitrogen, Oxygen, Air - 1.4
% Methane - 1.32
% Carbon dioxide - 1.28
% Propane - 1.13

ox_gamma = 1.667; % Specific heat ratio of OXIDIZER ullage gas
fu_gamma = 1.667; % Specific heat ratio of FUEL ullage gas

% initial properties
oxP0 = 300; % psia initial oxidizer tank pressure
oxV0 = 3.967; % ft3 initial oxidizer ullage volume
oxmdot0 = 11.464; % lbm/s initial oxidizer mass flow rate
fuP0 = 300; % psia initial fuel tank pressure
fuV0 = 2.894; % ft3 initial fuel ullage volume
fumdot0 = 6.034; % lbm/s initial fuel tank pressure
ofr0 = oxmdot0/fumdot0 % initial O/F ratio
oxrho = 68; % lbm/ft3 oxidizer density
furho = 49; % lbm/ft3 fuel density

T0 = 4435.744; % lbf, initial thrust
Pc0 = 120; % psia initial chamber pressure

% flow constants
oxbeta = Pc0/oxmdot0 % psia*s/lbm ox chamber pressure constant
oxalpha = oxmdot0/sqrt(oxP0-Pc0) % lbm/s/sqrt(psia) ox injector constant
fubeta = Pc0/fumdot0 % psia*s/lbm fuel chamber pressure constant
fualpha = fumdot0/sqrt(fuP0-Pc0) % lbm/s/sqrt(psia) fuel injector constant

% end conditions
oxVf = 8.245; % ft3 final oxidizer ullage volume
fuVf = 6.019; % ft3 final fuel ullage volume
tmax = 60; % maximum runtime
minofr = 1.7; % minimum o/f ratio
maxofr = 2.1; % maximum o/f ratio

dt = 0.003
n_points = tmax * 1/dt
t = 0;
n = 1;

% Simulation record vectors
tnum = zeros(1,n_points);
oxPnum = fuPnum = oxVnum = fuVnum = oxmdotnum = fumdotnum = Pcnum = oxdPnum = fudPnum = ofrnum = fTnum = tnum;

% Current frame variables
oxV = oxV0;
fuV = fuV0;
oxP = oxP0;
fuP = fuP0;
oxmdot = oxmdot0;
fumdot = fumdot0;
ofr = ofr0;
Pc = Pc0;

% Simulation loop
% Exit conditions:
% - Exceeds size of preallocated vectors
% - Runs out of oxidizer or fuel (triggers warning if runs out of fuel first)
% - Exceeds maximum runtime
% - Leaves O/F ratio bounds (triggers warning)
% - Flow reversal (triggers warning) (Chamber pressure greater than feed pressure)
while n < n_points && oxV < oxVf && fuV < fuVf && t < tmax && minofr < ofr < maxofr && oxP > Pc && fuP > Pc
  oxdV = oxalpha * sqrt(oxP - Pc) / oxrho * dt; % ft3 ox change in ullage volume
  fudV = fualpha * sqrt(fuP - Pc) / furho * dt; % ft3 fuel change in ullage volume
  oxmdot = oxdV / dt * oxrho; % lbm/s ox mass flow rate
  fumdot = fudV / dt * furho; % lbm/s fuel mass flow rate
  if isothermal
    oxP = oxP0 * oxV0 / oxV; % psia new oxidizer pressure
    fuP = fuP0 * fuV0 / fuV; % psia new fuel pressure
  else
    oxP = oxP0*(oxV/(oxV+oxdV))^ox_gamma; % psia new oxidizer pressure
    fuP = fuP0*(fuV/(fuV+fudV))^fu_gamma; % psia new fuel pressure
  end
  oxV += oxdV; % ft3 new oxidizer ullage volume
  fuV += fudV; % ft3 new fuel ullage volume
  Pc = (fumdot * fubeta + oxmdot * oxbeta) / 2; % psia chamber pressure
  % Storing data in record vectors
  tnum(n) = t; % time, seconds
  oxPnum(n) = oxP; % oxidizer tank pressure, psia
  fuPnum(n) = fuP; % fuel tank pressure, psia
  oxVnum(n) = oxV; % oxidizer tank ullage volume, ft3
  fuVnum(n) = fuV; % fuel tank ullage volume, ft3
  oxmdotnum(n) = oxmdot; % oxidizer mass flow rate, lbm/s
  fumdotnum(n) = fumdot; % fuel mass flow rate, lbm/s
  ofrnum(n) = oxmdot/fumdot; % oxidizer/fuel mass flow rate ratio
  Pcnum(n) = Pc; % chamber pressure, psia
  oxdPnum(n) = oxP - Pc; % oxidizer injector pressure drop, psi
  fudPnum(n) = fuP - Pc; % fuel injector pressure drop, psi
  fTnum(n) = T0/Pc0 * Pc; % thrust, lbf
  t += dt; % incrementing time
  n += 1; % index increment
  % Printing reason for exiting loop
  if n >= n_points
    disp("Exit! Exceeded preallocated vector size (n_points)!")
  end
  if oxV >= oxVf
    disp("Exit! Ran out of oxidizer!")
  end
  if fuV >= fuVf
    disp("Exit! WARNING: Ran out of fuel!")
  end
  if t >= tmax
    disp("Exit! Exceeded maximum runtime (tmax)!")
  end
  if ofr < minofr || ofr > maxofr
    disp("Exit! WARNING: Left O/F ratio bounds!")
  end
  if oxP < Pc || fuP < Pc
    disp("Exit! WARNING: Flow reversal!")
  endif

end

disp("Final time:");
disp(t);

% Plot tank and chamber pressures over time
delete(figure(1));
figure(1);
hold on;
plot(tnum(1:(n-1)),oxPnum(1:(n-1)),"blue");
plot(tnum(1:(n-1)),fuPnum(1:(n-1)),"red");
plot(tnum(1:(n-1)),Pcnum(1:(n-1)),"green");
title("pressures vs time");
xlabel("seconds");
ylabel("psia");

% Plot O/F ratio over time
delete(figure(2));
figure(2);
plot(tnum(1:(n-1)),ofrnum(1:(n-1)));
title("O/F ratio vs time");
ylabel("O/F (wt)");
xlabel("seconds");

% Plot tank ullage volumes over time
delete(figure(3));
figure(3);
hold on;
plot(tnum(1:(n-1)),oxVnum(1:(n-1)),"blue");
plot(tnum(1:(n-1)),fuVnum(1:(n-1)),"red");
title("ullage volumes vs time");
xlabel("seconds");
ylabel("cubic feet");

% Plot fuel & oxidizer mass flow rates over time
delete(figure(4));
figure(4);
hold on;
plot(tnum(1:(n-1)),oxmdotnum(1:(n-1)),"blue");
plot(tnum(1:(n-1)),fumdotnum(1:(n-1)),"red");
title("mass flow rates vs time");
ylabel("pounds/second");
xlabel("seconds");

% Plot thrust over time
delete(figure(5));
figure(5);
hold on;
plot(tnum(1:(n-1)),fTnum(1:(n-1)),"blue");
title("Thrust vs time");
ylabel("pounds");
xlabel("seconds");

% Perform polynomial regression on thrust data
Tfit = polyfit(tnum(1:n-1),fTnum(1:n-1),4);
Tpred = polyval(Tfit,tnum(1:n-1));
% Plot regression polynomial onto thrust graph
plot(tnum(1:n-1),Tpred,"green");
% Print polynomial coefficients to command window
disp("Plotted regression predicted thrust in green!");
disp("Thrust as function of time:");
disp("Read as polynomial coefficients in decreasing order");
disp("at^4 + bt^3 + ct^2 + dt + e");
disp(Tfit);


