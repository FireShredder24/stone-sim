% Pressure regulated rocket system simulator
% (c) John Nguyen 2025
% MIT License

% Common specific heat ratios:
% Helium - 1.667
% Hydrogen - 1.41
% Nitrogen, Oxygen, Air - 1.4
% Methane - 1.32
% Carbon dioxide - 1.28
% Propane - 1.13


% Real waterflow data
data = [
221.376243
223.2333
222.309411
222.539263
227.236833
225.066109
225.226223
225.884365
225.615293
225.677547
224.490279
222.695793
222.48568
221.331938
219.710427
217.048871
214.74708
208.981056
203.318903
200.827953
200.784832
197.494499
191.331052
186.966859
184.83417
181.996056
178.950835
174.59173
171.605467
167.613929
168.75299
163.000751
160.207232
157.880268
153.187207
154.000347
150.880011
147.955396
145.192108
144.51122
142.517402
140.671703
139.269265
137.359231
131.313066
130.759085
123.790548
116.765132
108.257669
101.984065
];
data = data + 14.7;
datatime = 0:0.5:(length(data)/2-0.5);

hp_data = [
395.629104
380.991364
366.803816
353.68294
343.123317
331.746456
317.167682
304.184085
291.953386
282.497666
270.837993
260.263766
249.460604
239.834916
230.89784
223.282906
214.579508
209.138799
201.793594
198.772662
196.371314
192.13125
188.55838
185.469915
181.465656
178.762407
174.912075
172.176655
169.476109
167.579664
165.233084
165.01634
161.287158
157.834578
155.7071
154.95748
151.430437
148.980847
147.458637
144.354675
142.353724
139.728503
137.998777
135.313698
131.368376
129.710999
123.738851
116.856287
110.706998
103.820869
];

hp_data = hp_data + 14.7
hp_datatime = 0:0.5:(length(hp_data)/2-0.5);



ox_gamma = 1.4; % Specific heat ratio of OXIDIZER ullage gas

% initial properties
resP0 = 415; % psia initial high pressure reservoir pressure
resV = 416/1728 % ft3 reservoir volume

oxP0 = 221.376+14.7; % psia initial oxidizer tank pressure
oxV0 = (736-677)/1728; % ft3 initial oxidizer ullage volume

oxrho = 62.34; % lbm/ft3 oxidizer density

% flow constants

oxCv = 0.7;

% end conditions
oxVf = 736/1728; % ft3 final oxidizer ullage volume
tmax = 60; % maximum runtime

enableplot = false

dt = 0.003
n_points = ceil(tmax * 1/dt)
t = 0;
n = 1;

% Simulation record vectors
tnum = zeros(1,n_points);
oxPnum = oxVnum = oxmdotnum = oxdPnum = resPnum = tnum;

% Current frame variables
oxV = oxV0;
oxP = oxP0;
resP = resP0;
mass = 0;
Pa = 14.7;

% Simulation loop
% Exit conditions:
% - Exceeds size of preallocated vectors
% - Runs out of oxidizer or fuel (triggers warning if runs out of fuel first)
% - Exceeds maximum runtime
% - Flow reversal (triggers warning) (Chamber pressure greater than feed pressure)
while n < n_points && oxV < oxVf && t < tmax && oxP > Pa
  oxdV = oxCv * sqrt((oxP - Pa)/(oxrho/62.34)) / 7.48 / 60 * dt; % ft3 ox change in ullage volume
  oxV2 = oxV + oxdV; % ft3 new oxidizer ullage volume

  % Liquid flow rates out are calculated from using the change in volume
  oxmdot = oxdV / dt * oxrho; % lbm/s ox mass flow rate
  mass += oxmdot * dt; % lbm, ox mass flowed

  % New tank pressures are found using blowdown methods first
  oxP2 = oxP * oxV / (oxV+oxdV); % psia new oxidizer pressure

  % Then the regulators kick in and keep pressures high, until the reservoir runs dry
  if resP > oxP0 && (resP-oxP0)*resV>(oxP2-oxP0)*oxV2
    % First case is when the reservoir has enough to keep constant pressure
    oxP2 = oxP0;
    % Derived from Boyle's Law relation:
    % ResP1V + OxP1V1 = ResP2V + OxP2V2
    % this distributes the gas among the reservoir and tank
    resP2 = (resP*resV + oxP*oxV - oxP0*oxV2) / resV;
  else
    % Second case is when the reservoir can't keep up
    % Derived from Boyle's Law relation:
    % ResP1V + OxP1V1 = ResP2V + OxP2V2
    % Which because reservoir and tank pressures are now equal simplifies to:
    % ResP1V + OxP1V1 = ResP2*(ResV+OxV2) = OxP2*(ResV+OxV2)
    oxP2 = resP2 = (resP*resV + oxP*oxV)/(resV+oxV2);
  endif

  oxP = oxP2;
  resP = resP2;
  oxV += oxdV; % ft3 new oxidizer ullage volume
  % Storing data in record vectors
  tnum(n) = t; % time, seconds
  oxPnum(n) = oxP; % oxidizer tank pressure, psia
  oxVnum(n) = oxV; % oxidizer tank ullage volume, ft3
  oxmdotnum(n) = oxmdot; % oxidizer mass flow rate, lbm/s
  resPnum(n) = resP;
  oxdPnum(n) = oxP - Pa; % oxidizer injector pressure drop, psi
  t += dt; % incrementing time
  n += 1; % index increment
  % Printing reason for exiting loop
  if n >= n_points
    disp("Exit! Exceeded preallocated vector size (n_points)!")
  end
  if oxV >= oxVf
    disp("Exit! Ran out of oxidizer!")
  end
  if t >= tmax
    disp("Exit! Exceeded maximum runtime (tmax)!")
  end
  if oxP < Pa
    disp("Exit! WARNING: Flow reversal!")
  end

end

disp("Final time: (s)");
disp(t);
disp("Total mass flowed: (lbm)");
disp(mass);
disp("Average mass flow rate: (lbm/s)");
disp(mass/t);

% Plot tank and chamber pressures over time
delete(figure(1));
figure(1);
hold on;
plot(tnum(1:(n-1)),oxPnum(1:(n-1)),"blue");
plot(tnum(1:(n-1)),resPnum(1:(n-1)),"green");
plot(datatime,data,"o","linestyle","none");
plot(hp_datatime,hp_data,"+","linestyle","none");
title("Pressures vs Time");
legend("Sim Tank","Sim High Press","Real Tank Data","Real HP Data");
xlabel("Time - Seconds");
ylabel("Pressure - psia");

if enableplot
% Plot tank ullage volumes over time
delete(figure(3));
figure(3);
hold on;
plot(tnum(1:(n-1)),oxVnum(1:(n-1)),"blue");
title("ullage volumes vs time");
xlabel("seconds");
ylabel("cubic feet");

% Plot fuel & oxidizer mass flow rates over time
delete(figure(4));
figure(4);
hold on;
plot(tnum(1:(n-1)),oxmdotnum(1:(n-1)),"blue");
title("mass flow rates vs time");
ylabel("pounds/second");
xlabel("seconds");
end
