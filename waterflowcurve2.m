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
% TANK PRESSURE DATA
data = [
226.148032
225.359718
224.616288
224.221525
223.974938
223.646619
223.675925
222.964547
222.474004
222.471027
221.297834
221.170321
220.197964
218.829022
214.920373
208.119709
202.093312
197.235629
193.300044
188.844292
185.074485
181.359475
177.689031
174.181417
171.112844
167.738118
165.378553
162.750234
159.5036
156.928546
154.485485
150.844405
150.201234
147.578898
145.64372
143.062684
141.707239
139.62689
137.607032
136.153698
134.028176
132.467416
130.488828
129.309075
127.287092
];
data = data + 14.7;
datatime = 0:0.5:(length(data)/2-0.5);

% HIGH PRESS DATA
hp_data = [380.942547
366.534087
352.34145
338.647477
325.616817
312.487029
299.969925
288.243265
276.064661
264.856301
253.284459
242.80016
232.416169
222.851804
215.650089
210.441275
205.804379
200.724192
196.509744
191.871353
187.78671
183.989362
180.379291
176.956178
173.563426
170.708034
167.669564
164.81538
162.383068
159.50804
157.008426
154.603772
152.733834
149.958251
147.975183
145.509493
143.779796
141.818434
140.165843
138.499826
136.47924
134.685804
133.054029
131.757791
130.078679
];

hp_data = hp_data + 14.7;
hp_datatime = 0:0.5:(length(hp_data)/2-0.5);

% INLET PRESSURE DATA
inlet_data = [216.049044
217.191909
215.515008
214.111384
213.031111
211.724037
210.398488
210.907199
209.411197
209.933289
208.620575
208.814569
206.907615
206.352715
201.728819
195.618763
190.430441
185.484978
181.864777
177.932246
174.353168
171.084002
167.167875
164.151438
160.680538
158.247002
155.212665
153.61476
150.769956
147.641458
146.038775
143.56918
139.038266
139.428929
136.765053
134.93675
133.462821
131.618748
129.676599
128.347481
126.353215
124.801786
123.280143
126.761202
126.450974
];

inlet_data = inlet_data + 14.7;
inlet_datatime = 0:0.5:(length(inlet_data)/2-0.5);

ox_gamma = 1.4; % Specific heat ratio of OXIDIZER ullage gas

% initial properties
resP0 = hp_data(1); % psia initial high pressure reservoir pressure
resV = 420/1728; % ft3 reservoir volume

oxP0 = data(1); % psia initial oxidizer tank pressure
oxV0 = (736-677-40.5)/1728; % ft3 initial oxidizer ullage volume

oxrho = 62.34; % lbm/ft3 oxidizer density

tankPT_thru_diameter = 0.75; % in
tankPT_flow_area = tankPT_thru_diameter^2/4*pi/144; % ft2

highpressPT_thru_diameter = 0.305; % in
highpressPT_flow_area = highpressPT_thru_diameter^2/4*pi/144; % ft2

inletPT_thru_diameter = 0.43; % in
inletPT_flow_area = inletPT_thru_diameter^2/4*pi/144; % ft2

% flow constants

totalCv = 0.73; % gal/min/sqrt(psi) Overall flow coefficient

CvRatio = 4; % Ratio between lower fluids and orifice Cv's

oxOrificeCv = sqrt(totalCv^2*(1+1/CvRatio^2)) % gal/min/sqrt(psi) Orifice flow coefficient
oxLowerCv = oxOrificeCv*CvRatio % gal/min/sqrt(psi) Lower fluids flow coefficient

totalCv = sqrt(1/(1/oxLowerCv^2+1/oxOrificeCv^2)); % gal/min/sqrt(psi) Overall flow coefficient

% end conditions
oxVf = (736)/1728; % ft3 final oxidizer ullage volume
tmax = 60; % maximum runtime

enableplot = false;

dt = 0.005;
n_points = ceil(tmax * 1/dt);
t = 0;
n = 1;

% Simulation record vectors
tnum = zeros(1,n_points);
oxPnum = tnum;
oxVnum = tnum;
oxmdotnum = tnum;
oxdPnum = tnum;
resPnum = tnum;
inletPnum = tnum;

% Current frame variables
oxV = oxV0;
oxP = oxP0;
Q = 25.88/22/62.34; % ft3/s, initial flowrate estimate
inletP = 0;
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
  tank_gas_density = oxP/10.731*28/530; % lbm/ft3 Density of gas in the propellant tank (ideal gas law)
  tank_dynamic_pressure = tank_gas_density*(Q/tankPT_flow_area)^2/2/144/32.174; % psid, dynamic pressure at tank PT
  tank_P_eff = oxP + tank_dynamic_pressure; % psia, effective tank pressure after gas comes to "rest" in the tank
  Q = totalCv * sqrt((tank_P_eff - Pa)/(oxrho/62.34)) / 7.48052 / 60; % ft3/s, liquid flowrate through downstream constrictions

  oxdV = Q * dt; % ft3 ox change in ullage volume
  oxV2 = oxV + oxdV; % ft3 new tank ullage volume

  inlet_dynamic_pressure = oxrho*(Q/inletPT_flow_area)^2/2/144/32.174; % psi, dynamic pressure at inlet PT tee
  inletP = ((Q*7.48052*60)^2/oxOrificeCv^2)*(oxrho/62.34) + Pa; % psia, simulated inlet PT reading

  % Liquid flow rates out are calculated from using the change in volume
  oxmdot = oxdV / dt * oxrho; % lbm/s, ox mass flow rate
  mass = mass + oxmdot * dt; % lbm, ox mass flowed

  % New tank pressures are found using blowdown methods first
  oxP2 = tank_P_eff * oxV / (oxV+oxdV); % psia, new tank pressure
  
  res_gas_density = resP/10.731*28/530; % lbm/ft3, Density of the gas in the reservoir (approximate)
  res_dynamic_pressure = res_gas_density*(Q/tankPT_flow_area)^2/2/144/32.174; % psi, dynamic pressure at the tee for the high press PT
  resP_eff = resP + res_dynamic_pressure; % psia, static pressure in the reservoir

  % Then the regulators kick in and keep pressures high, until the reservoir runs dry
  if resP_eff > oxP0 && (resP_eff-oxP0)*resV>(oxP2-oxP0)*oxV2
    % First case is when the reservoir has enough to keep constant pressure
    oxP2 = oxP0;
    % Derived from Boyle's Law relation:
    % ResP1V + OxP1V1 = ResP2V + OxP2V2
    % this distributes the gas among the reservoir and tank
    resP2 = (resP_eff*resV + oxP*oxV - oxP0*oxV2) / resV;
  else
    % Second case is when the reservoir can't keep up
    % Derived from Boyle's Law relation:
    % ResP1V + OxP1V1 = ResP2V + OxP2V2
    % Which because reservoir and tank pressures are now equal simplifies to:
    % ResP1V + OxP1V1 = ResP2*(ResV+OxV2) = OxP2*(ResV+OxV2)
    oxP2 = (resP_eff*resV + oxP*oxV)/(resV+oxV2);
    resP2 = oxP2;
  end

  oxP = oxP2 - tank_dynamic_pressure; % psia, new simulated tank PT reading
  resP = resP2 - res_dynamic_pressure; % psia, new simulated high press PT reading
  oxV = oxV + oxdV; % ft3 new oxidizer ullage volume
  % Storing data in record vectors
  tnum(n) = t; % time, seconds
  oxPnum(n) = oxP; % oxidizer tank pressure, psia
  oxVnum(n) = oxV; % oxidizer tank ullage volume, ft3
  oxmdotnum(n) = oxmdot; % oxidizer mass flow rate, lbm/s
  resPnum(n) = resP;
  inletPnum(n) = inletP;
  oxdPnum(n) = oxP - Pa; % oxidizer injector pressure drop, psi
  t = t + dt; % incrementing time
  n = n + 1; % index increment
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
plot(tnum(1:(n-1)),inletPnum(1:(n-1)),"red");
plot(datatime,data,"o","linestyle","none");
plot(hp_datatime,hp_data,"+","linestyle","none");
plot(inlet_datatime,inlet_data,"x","linestyle","none");
title("Pressures vs Time");
legend("Sim Tank","Sim High Press","Sim Inlet","Real Tank Data","Real High Press Data", "Real Inlet Data");
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
