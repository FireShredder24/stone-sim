# By John Nguyen <knightofthealtar64@gmail.com>
# Python 3.10.11 with VPython 7.6.4
# MIT License

import math
import random
from vpython import *


# UTILITY FUNCTIONS

# yaw & pitch angle to 3d xyz unit vector
def yp_unit(yaw, pitch):
    return vec(
        cos(yaw) * cos(pitch),  # x
        sin(pitch),  # y
        sin(yaw) * cos(pitch)  # z
    )


# ENVIRONMENT FUNCTIONS
# bulk modulus of elasticity of air (used for speed of sound)
bulk_mod = 1.01325 * 10 ** 5  # Pa

# Specific gas constant of air (R/M)
r_spec = 287.05  # J/(kg * K)


# air density function, based on International Standard Atmosphere
def rho(y):  # kg/m^3 air pressure
    dens = 1.23
    if y <= 11000:
        dens = -(1.2985 - 0.3639) / 11610 * (y + 610) + 1.2985
    elif 11000 < y <= 20000:
        dens = -(0.3639 - 0.088) / (20000 - 11000) * (y - 11000) + 0.3639
    elif 20000 < y <= 32000:
        dens = -(0.088 - 0.0132) / (32000 - 20000) * (y - 20000) + 0.088
    elif 32000 < y <= 47000:
        dens = -(0.0132 - 0.002) / (47000 - 32000) * (y - 32000) + 0.0132
    else:
        dens = 0
    return dens


# air temperature function, from ISA
def temp(y):  # K temperature
    T = 273
    if y <= 11000:
        T = 273 + 19 - 6.5 * y / 1000  # -6.5 K per km
    elif 11000 < y <= 20000:
        T = 273 - 56.5 - 0 * y / 1000  # +0.0 K per km
    elif 20000 < y <= 32000:
        T = 273 - 56.5 + 1.0 * y / 1000  # +1.0 K per km
    elif 32000 < y <= 47000:
        T = 273 - 44.5 + 2.8 * y / 1000  # +2.8 K per km
    else:
        T = 273 - 44.5
    return T


# static pressure function
def P(y):
    return rho(y) * r_spec * temp(y)


# speed of sound function
def c(y):
    return sqrt(bulk_mod / rho(y))


g_0 = -9.0866  # m/s^2, gravitational acceleration at MSL


# Gravitational acceleration function
def g(y):
    return g_0 * (6371009 / (6371009 + y)) ** 2


# wind profile class declaration
class WindProfile:
    def __init__(self, name, mu, sigma, angle_mu, angle_sigma, step):
        self.name = name

        # GENERATION PARAMETERS
        self.mu = mu  # m/s, mean windspeed
        self.sigma = sigma  # m/s, std. dev. of windspeed
        self.angle_mu = angle_mu  # rad, mean change in angle per band
        self.angle_sigma = angle_sigma  # rad, std. dev. of change in angle per band
        self.step = step  # altitude band size

        alt = -50
        self.altSet = list()
        while alt < 47000:
            self.altSet.append(alt)
            alt += step
        del alt

        self.speedSet = list()
        self.angleSet = list()
        for y in self.altSet:
            speed = random() * self.sigma + self.mu
            angle = random() * self.angle_sigma + self.angle_mu
            self.speedSet.append(speed)
            self.angleSet.append(angle)

        print("speed set:", self.speedSet)
        print("angle set:", self.angleSet)
        # end constructor

    # TODO: Interpolate between wind bands
    def wind(self, altitude):
        index = int((altitude + 50) / self.step)
        return yp_unit(self.angleSet[index], 0) * self.speedSet[index]


# Physics object class declaration
class FreeRocket:
    # TODO: Include switches for all graphs
    # constructor
    def __init__(self, name, pos, yaw, pitch, roll, v_0, ymi, pmi, rmi, cp, cd, A, cd_s, A_s, chute_cd, chute_A,
                 drogue_cd,
                 drogue_A, cg, dry_mass, fuel_mass, thrust, t0, wind):
        self.name = name
        # FUNDAMENTAL VECTORS
        self.pos = pos  # m, 3D cartesian position
        self.rot = vec(yaw, pitch, roll)  # rad, rotational position
        self.v_0 = v_0  # m/s, initial velocity magnitude
        self.I_0 = vec(ymi, pmi, rmi)  # kg*m^2, Mass moments of inertia
        # AERODYNAMIC PROPERTIES
        # TODO: Make cp and cg 3D
        self.cp = cp  # m, center of pressure position aft of nose cone
        self.cd = cd  # frontal drag coefficient
        self.A = A  # m^2, frontal reference area
        self.cd_s = cd_s  # side drag coefficient (airflow parallel to yaw/pitch axis)
        self.A_s = A_s  # side reference area " "
        self.chute_cd = chute_cd  # parachute drag coefficient
        self.chute_A = chute_A  # m^2, parachute reference area
        self.drogue_cd = drogue_cd  # drogue chute cd
        self.drogue_A = drogue_A  # m^2, drogue chute ref. area
        # MASS PROPERTIES
        self.cg = cg  # m, center of mass position aft of nose cone
        self.dry_mass = dry_mass  # kg, total dry mass
        self.fuel_mass = fuel_mass * 0.453  # lb to kg
        self.mass = self.dry_mass + self.fuel_mass  # kg, total initial mass
        # PROPULSION PROPERTIES
        self.thrust = thrust  # thrust function of time
        self.wind = wind  # Wind profile

        t = 0  # s, thrust time variable
        dt = 0.01  # s, thrust time step
        self.J = 0  # total impulse
        while self.thrust(t) > 0:  # integration of thrust
            self.J += self.thrust(t)
            t += dt
        self.bt = t  # s, motor burn time
        self.t0 = t0  # s, ignition time
        self.t1 = self.t0 + self.bt  # s, burnout time
        # MOMENTUM PROPERTIES
        self.p = self.mass * self.v_0 * yp_unit(self.rot.x, self.rot.y)  # linear momentum
        self.v = self.p / self.mass  # Linear velocity
        self.drot = vec(0, 0, 0)  # rad/s, yaw/pitch/roll rate
        self.L = vec(  # Angular momentum, yaw/pitch/roll
            self.I_0.x * self.drot.x,
            self.I_0.y * self.drot.y,
            self.I_0.z * self.drot.z
        )
        self.drogue = False
        self.main_chute = False

        # TODO: More graphs and graphing options

        self.flight_side_graph = graph(title="Side Profile of " + self.name, xtitle="Downrange Distance",
                                       ytitle="Altitude (m)")
        self.flight_side = gcurve(graph=self.flight_side_graph, color=color.green)

        self.alpha_graph = graph(title="AoA of " + self.name, xtitle="t", ytitle="alpha, radians")
        self.alpha_graph_curve = gcurve(graph=self.alpha_graph, color=color.red)

        self.moment_graph = graph(title="Moments on " + self.name, xtitle="t", ytitle="Nm")
        self.moment_graph_y = gcurve(graph=self.moment_graph, color=color.red)
        self.moment_graph_p = gcurve(graph=self.moment_graph, color=color.green)
        self.moment_graph_r = gcurve(graph=self.moment_graph, color=color.blue)

        self.velocity_graph = graph(title="Velocity of " + self.name, xtitle="t", ytitle="m/s")
        self.velocity_graph_x = gcurve(graph=self.velocity_graph, color=color.red)
        self.velocity_graph_y = gcurve(graph=self.velocity_graph, color=color.green)
        self.velocity_graph_z = gcurve(graph=self.velocity_graph, color=color.blue)
        self.velocity_graph_total = gcurve(graph=self.velocity_graph, color=color.black)

        self.v_max = 0  # m/s, maximum absolute velocity
        self.y_max = 0  # m, maximum altitude
        self.a_max = 0  # m/s^2, maximum acceleration
        self.g_max = 0  # G, maximum acceleration
        # end of constructor

    # Drag coefficient estimator
    def cd_alpha(self, alpha):  # Generates parabola which interpolates b/w cd and cd_s from -pi/2 to +pi/2
        return alpha ** 2 / ((pi / 2) ** 2 / (self.cd_s - self.cd)) + self.cd

    # Reference area estimator
    def A_alpha(self, alpha):  # Generates parabola which interpolates b/w A and A_s from -pi/2 to +pi/2
        return alpha ** 2 / ((pi / 2) ** 2 / (self.A_s - self.A)) + self.A

    # Simulation 
    def simulate(self, t, dt):
        # FORCE & MOMENT COMPONENTS

        # BODY DRAG
        heading = yp_unit(self.rot.x, self.rot.y)  # 3d linear unit vector of vehicle orientation
        airflow = self.v + self.wind.wind(self.pos.y)  # 3d linear vector of oncoming airstream (reversed)
        alpha = math.acos(heading.dot(airflow.hat))  # rad, angle of attack
        self.alpha_graph_curve.plot(t, alpha)
        cd = self.cd_alpha(alpha)  # drag coefficient
        A = self.A_alpha(alpha)  # m^2, reference area
        f_drag = self.v.mag ** 2 * rho(self.pos.y) * cd * A / 2 * -self.v.hat  # N, body drag force

        M_drag = f_drag.cross(vec(0, self.cg, 0) - vec(0, self.cp, 0))  # Torque generated by drag force vector

        # DROGUE PARACHUTE DRAG
        if self.drogue:
            f_drogue = self.v.mag ** 2 * rho(self.pos.y) * self.drogue_cd * self.drogue_A / 2 * -self.v.hat
        else:
            f_drogue = vec(0, 0, 0)

        # MAIN PARACHUTE DRAG
        if self.main_chute:
            f_chute = self.v.mag ** 2 * rho(self.pos.y) * self.chute_cd * self.chute_A / 2 * -self.v.hat
        else:
            f_chute = vec(0, 0, 0)

        # GRAVITY FORCE
        f_grav = vec(0, self.mass * g(self.pos.y), 0)

        # THRUST VECTOR
        f_thrust = vec(0, 0, 0)
        if self.t0 <= t <= self.t1:
            f_thrust = self.thrust(t - self.t0) * yp_unit(self.rot.x, self.rot.y)

        # TODO: Update mass based on consumed fuel

        # TOTAL NET FORCE
        f_net = f_grav + f_thrust + f_drag + f_chute + f_drogue

        # TOTAL NET MOMENT
        M_net = M_drag

        self.p = self.p + f_net * dt  # incrementing linear momentum
        self.v = self.p / self.mass  # calculating velocity
        self.pos = self.pos + self.v * dt  # incrementing position

        self.L += M_net * dt  # incrementing angular momentum
        self.drot = self.L / self.mass  # calculating rotation rate
        self.rot = self.rot + self.drot * dt  # incrementing orientation

        if self.v.y < 5 and not self.drogue:
            self.drogue = True
            print("Drogue deployment at T+", "{:3.2f}".format(t), " at Altitude: ", "{:3.0f}".format(self.pos.y), " & Speed: ", "{:3.2f}".format(self.v.mag))

        if self.pos.y < 150 and self.drogue and not self.main_chute:
            self.main_chute = True
            print("Main chute deployment at T+", "{:3.2f}".format(t), " at Altitude: ", "{:3.0f}".format(self.pos.y), " & Speed: ", "{:3.2f}".format(self.v.mag))

        self.flight_side.plot(payload.pos.x, payload.pos.y)
        self.moment_graph_y.plot(t, M_net.x)
        self.moment_graph_p.plot(t, M_net.y)
        self.moment_graph_r.plot(t, M_net.z)
        self.velocity_graph_x.plot(t, self.v.x)
        self.velocity_graph_y.plot(t, self.v.y)
        self.velocity_graph_z.plot(t, self.v.z)
        self.velocity_graph_total.plot(t, self.v.mag)
        
        if self.pos.y > self.y_max:
            self.y_max = self.pos.y
        if self.v.mag > self.v_max:
            self.v_max = self.v.mag
        a = f_net.mag / self.mass
        if a > self.a_max:
            self.a_max = a
        if a / g_0 > self.g_max:
            self.g_max = a / g_0
        # end of simulation method
    # TODO: Implement after-flight report method
    # end of class def


# payload/e-bay thrust function
def payload_thrust(t):
    return -77 * (t - 0.2) ** 2 + 300


wind_1 = WindProfile("FAR", 5, 3, pi / 4, pi, 500)

# Alpha Phoenix on I-300
payload = FreeRocket(
    "Alpha Phoenix",
    vec(0, 1, 0),  # 3D linear position
    0,  # yaw
    88 * pi / 180,  # pitch
    0,  # roll
    5,  # v_0
    0.0715,  # yaw moment of inertia
    0.0715,  # pitch moment of inertia
    0.0012,  # roll moment of inertia
    0.152,  # center of pressure position aft of nose
    0.6,  # frontal drag coefficient
    0.015,  # frontal reference area
    1.5,  # side drag coefficient (orthogonal to long. axis)
    0.05,  # side reference area (orthogonal to long. axis)
    0.8,  # main parachute drag coefficient
    0.4,  # main parachute ref. area
    0.8,  # drogue chute drag coefficient
    0.1,  # drogue chute ref. area
    0.142,  # cg position aft of nose
    1.1,  # dry mass
    0.220,  # fuel mass
    payload_thrust,  # thrust function
    0,  # ignition time
    wind_1
)
print("payload p:", payload.p)
print("payload position: ", payload.pos, " rotation: ", payload.rot, " moments of inertia (YPR): ", payload.I_0)

time = 0
dtime = 0.05

print("----BEGIN SIMULATION----")
while time < 120 and payload.pos.y > 0:
    payload.simulate(time, dtime)

    time += dtime
print("----END SIMULATION----")

# TODO call summary stat method of each FreeRocket
