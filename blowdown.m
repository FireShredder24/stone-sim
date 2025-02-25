% Blowdown rocket system simulator
% (c) John Nguyen 2025
% MIT License

% initial properties
oxP0 = 300 % psia initial oxidizer tank pressure
oxV0 = 3.967 % ft3 initial oxidizer ullage volume
oxmdot0 = 11.464 % lbm/s initial oxidizer mass flow rate
fuP0 = 300 % psia initial fuel tank pressure
fuV0 = 2.894 % ft3 initial fuel ullage volume
fumdot0 = 6.034 % lbm/s initial fuel tank pressure
ofr0 = oxmdot0/fumdot0 % initial O/F ratio
oxrho = 68 % lbm/ft3 oxidizer density
furho = 49 % lbm/ft3 fuel density

T0 = 4435.744 % lbf, initial thrust
Pc0 = 120 % psia initial chamber pressure

% flow constants
oxbeta = Pc0/oxmdot0 % psia*s/lbm ox chamber pressure constant
oxalpha = oxmdot0/sqrt(oxP0-Pc0) % lbm/s/sqrt(psia) ox injector constant
fubeta = Pc0/fumdot0 % psia*s/lbm fuel chamber pressure constant
fualpha = fumdot0/sqrt(fuP0-Pc0) % lbm/s/sqrt(psia) fuel injector constant

% end conditions
oxVf = 8.245 % ft3 final oxidizer ullage volume
fuVf = 6.019 % ft3 final fuel ullage volume
tmax = 60 % maximum runtime
minofr = 1.7 % minimum o/f ratio
maxofr = 2.1 % maximum o/f ratio

dt = 0.003
n_points = tmax * 1/dt
t = 0
n = 1

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

while n < n_points && oxV < oxVf && fuV < fuVf && t < tmax && minofr < ofr < maxofr && oxP > Pc && fuP > Pc
  oxdV = oxalpha * sqrt(oxP - Pc) / oxrho * dt; % ft3 ox change in ullage volume
  fudV = fualpha * sqrt(fuP - Pc) / furho * dt; % ft3 fuel change in ullage volume
  oxmdot = oxdV / dt * oxrho; % lbm/s ox mass flow rate
  fumdot = fudV / dt * furho; % lbm/s fuel mass flow rate
  oxV += oxdV; % ft3 new oxidizer ullage volume
  fuV += fudV; % ft3 new fuel ullage volume
  oxP = oxP0 * oxV0 / oxV; % psia new oxidizer pressure
  fuP = fuP0 * fuV0 / fuV; % psia new fuel pressure
  Pc = (fumdot * fubeta + oxmdot * oxbeta) / 2; % psia chamber pressure
  t += dt; % incrementing time
  % Storing data in record vectors
  tnum(n) = t;
  oxPnum(n) = oxP;
  fuPnum(n) = fuP;
  oxVnum(n) = oxV;
  fuVnum(n) = fuV;
  oxmdotnum(n) = oxmdot;
  fumdotnum(n) = fumdot;
  ofrnum(n) = oxmdot/fumdot;
  Pcnum(n) = Pc;
  oxdPnum(n) = oxP - Pc;
  fudPnum(n) = fuP - Pc;
  fTnum(n) = T0/Pc0 * Pc; % lbf, thrust
  n += 1; % Incrementing index vector
  % Printing reason for exiting loop
  if n >= n_points
    disp("exceeded n_points")
  end
  if oxV >= oxVf
    disp("ran out of oxidizer")
  end
  if fuV >= fuVf
    disp("ran out of fuel")
  end
  if t >= tmax
    disp("exceeded tmax")
  end
  if ofr < minofr || ofr > maxofr
    disp("left O/F ratio bounds")
  end
  if oxP < Pc || fuP < Pc
    disp("flow reversal!")
  endif

end

disp("Final time:");
disp(t);

delete(figure(1));
figure(1);
hold on;
plot(tnum,oxPnum,"blue");
plot(tnum,fuPnum,"red");
plot(tnum,Pcnum,"green");
title("pressures vs time");
xlabel("seconds");
ylabel("psia");

delete(figure(2));
figure(2);
plot(tnum,ofrnum);
title("O/F ratio vs time");
ylabel("O/F (wt)");
xlabel("seconds");

delete(figure(3));
figure(3);
hold on;
plot(tnum,oxVnum,"blue");
plot(tnum,fuVnum,"red");
title("ullage volumes vs time");
xlabel("seconds");
ylabel("cubic feet");

delete(figure(4));
figure(4);
hold on;
plot(tnum,oxmdotnum,"blue");
plot(tnum,fumdotnum,"red");
title("mass flow rates vs time");
ylabel("pounds/second");
xlabel("seconds");

delete(figure(5));
figure(5);
hold on;
plot(tnum,fTnum,"blue");
title("Thrust vs time");
ylabel("pounds");
xlabel("seconds");
Tfit = polyfit(tnum(1:n-1),fTnum(1:n-1),4);
Tpred = polyval(Tfit,tnum(1:n-1));
plot(tnum(1:n-1),Tpred,"green");
disp("Plotted regression predicted thrust in green!");
disp("Thrust as function of time");
disp(Tfit);


