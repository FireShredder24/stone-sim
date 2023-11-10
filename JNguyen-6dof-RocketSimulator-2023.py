# By John Nguyen <knightofthealtar64@gmail.com>
# Python 3.10.11 with VPython 7.6.4
# MIT License

import math
import random
import numpy as np
from vpython import *


# UTILITY FUNCTIONS

# yaw & pitch angle to 3d xyz unit vector
def yp_unit(yaw, pitch):
    return vec(
        cos(yaw) * cos(pitch),  # x
        sin(pitch),  # y
        sin(yaw) * cos(pitch)  # z
    )

# Common number string formatting style
def d2(num):
    return "{:.2f}".format(num)

# ENVIRONMENT FUNCTIONS
# bulk modulus of elasticity of air (used for speed of sound)
bulk_mod = 1.01325 * 10 ** 5  # Pa

# Specific gas constant of air (R/M)
r_spec = 287.05  # J/(kg * K)


# air density function, based on International Standard Atmosphere
def rho(y):  # kg/m^3 air pressure
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


G = 6.674e-11  # Nm^2/kg^2, Universal gravitational constant

g_0 = 9.80665 # m/s^2, Earth standard gravity

M_E = 6e24 # kg, Earth mass

# Gravitational force function
# Takes altitude from mean sea level
def g(y: float, mass: float):
    return G * M_E * mass / (y+6371000)**2


# wind profile class declaration
class WindProfile:
    def __init__(self, name: str, mu: float, sigma: float, angle_mu: float, angle_sigma: float, step: float, print_debug: bool):
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

        if print_debug:
            print(f"speed set: {self.speedSet}")
            print(f"angle set: {self.angleSet}")
        # end constructor

    def wind(self, altitude):
        index = (altitude + 50) / self.step  # normalize to band number
        band_pos = index - floor(index)  # relative position within altitude band
        index = abs(int(index))  # convert to positive whole number for table index
        angle = self.angleSet[index]  # grab angle from wind band
        speed = self.speedSet[index]  # grab speed from wind band
        if band_pos <= 0.2:  # if in the lower 20% of wind band, interpolate to the avg. b/w current and previous bands
            angle_1 = self.angleSet[index - 1]  # angle from previous wind band
            speed_1 = self.speedSet[index - 1]  # speed from previous wind band
            # avg. of 2 wind bands (y-intercept) + difference b/w current band and avg. point (slope) * rel. pos (x)
            angle = (angle_1 + angle) / 2 + (angle - (angle_1 + angle) / 2) * band_pos / 0.2
            speed = (speed + speed_1) / 2 + (speed - (speed_1 + speed) / 2) * band_pos / 0.2
        if band_pos >= 0.8:  # if in the upper 20% of wind band, interpolate to the avg. b/w current and next bands
            angle_1 = self.angleSet[index + 1]  # angle from next wind band
            speed_1 = self.speedSet[index + 1]  # speed from next wind band
            # current wind band (y-intercept) + difference b/w current band and avg. point (slope) * rel. pos (x)
            angle = angle - (angle - (angle_1 + angle) / 2) * (band_pos - 0.8) / 0.2
            speed = speed - (speed - (speed_1 + speed) / 2) * (band_pos - 0.8) / 0.2

        return yp_unit(angle, 0) * speed  # unit vector of angle, times speed

# Fin set class declaration
# Determines lift properties of fins only.  Drag is handled by the overall rocket cd, A, cd_s, and A_s
class FinSet:
    def __init__(self, num_fins: int, center: vector, pos: vector, planform: float, stall_angle: float, ac_span: float, cl_pass):
        self.num_fins = num_fins # number of fins.  3 or 4 only supported at this time.
        self.fin_rad_pos = []
        for a in arange(0, 2*pi, 2*pi/num_fins):
            self.fin_rad_pos.append(a)
        print(f"Fin radial positions: {self.fin_rad_pos}")
        self.center = center # position of center of mass relative to parent rocket nose tip
        self.pos = pos # position of center of lift relative to parent rocket nose tip
        self.planform = planform # planform wing area of each individual fin
        self.stall_angle = stall_angle # maximum angle of attack before wing stall
        self.ac_span = ac_span # radial offset from rocket centerline to center of lift of each fin
        self.cl = cl_pass
        self.aoa_graph = graph(fast=False, title="Fin AoA", xtitle="t", ytitle="rad")
        self.aoa_curves = []
        for i in self.fin_rad_pos:
            self.aoa_curves.append(gcurve(graph=self.aoa_graph, label=f"Fin {i} AoA", color=color.blue))
        self.curve = gcurve(graph=self.aoa_graph, label="fin 0 aoa", color=color.red)

    def ac_pos(self, rot: vector, fin_index: int):
        if fin_index > self.num_fins - 1:
            print(f"Expected fin index b/w 0 and {self.num_fins-1} got {fin_index}")
            exit(1)
        return yp_unit(self.fin_rad_pos[fin_index] + rot.x, 0) * self.ac_span


    def lift(self, aoa: float, altitude: float, airspeed: float):
        return self.cl(aoa) * rho(altitude) * airspeed**2 / 2 * self.planform

    # Returns lift vector in rocket-centric coordinates.
    def lift_vec(self, rot: vector, fin_index: int, aoa: float, altitude: float, airflow: vector, roll: vector):
        L = self.lift(aoa, altitude, airflow.mag) * self.ac_pos(rot, fin_index).cross(roll).hat
        if mag(L + self.flow_perp_fin(rot, airflow, self.ac_pos(rot, fin_index), roll)) < mag(L):
            L *= -1
        return L

    # Returns the component of the airflow vector orthogonal to the fin plane.
    def flow_perp_fin(self, rot: vector, airflow: vector, ac: vector, roll:vector):
        f_II_ac = airflow.dot(ac.hat) * ac.hat # airflow vector projected onto the aerodynamic center pos vector
        f_II_roll = airflow.dot(roll.hat) * roll.hat # airflow vector projected onto the roll axis
        f_I_accg = airflow - (f_II_roll + f_II_ac) # remaining component of airflow vector (orthogonal to both ac pos & roll axis)
        return f_I_accg
 
    # Returns AoA of each fin.
    def aoa(self, rot: vector, fin_index: int, airflow: vector, roll: vector):
        ac = self.ac_pos(rot, fin_index) # aerodynamic center position vector
        f_II_ac = airflow.dot(ac.hat) * ac.hat # airflow vector projected onto the aerodynamic center pos vector
        f_II_roll = airflow.dot(roll.hat) * roll.hat # airflow vector projected onto the roll axis
        f_I_accg = airflow - (f_II_roll + f_II_ac) # remaining component of airflow vector (orthogonal to both ac pos & roll axis)
        alpha = atan(f_I_accg.mag/f_II_roll.mag) # right triangle with orthogonal component as opposite and roll component as adjacent
        return alpha

    # Iterates through each fin and plots its angle of attack
    def aoa_plot(self, t: float, rot:vector, airflow:vector, roll:vector):
        for idx, i in enumerate(self.aoa_curves):
            i.plot(t, self.aoa(rot, idx, airflow, roll))

    def total_lift_vec(self, rot: vector, airflow: vector, roll: vector, altitude: float):
        L_total = vec(0,0,0)
        for idx, i in enumerate(self.fin_rad_pos):
            aoa = self.aoa(rot, idx, airflow, roll)
            L = self.lift_vec(rot, idx, aoa, altitude, airflow, roll)
            L_total = L_total + L
        return L_total



# Physics object class declaration
class FreeRocket:
    # Graph switches
    position_graph_enable = False  
    rotation_graph_enable = False  
    velocity_graph_enable = True  
    rotation_rate_graph_enable = False  
    acceleration_graph_enable = True  
    moment_graph_enable = False  
    side_profile_enable = True  
    top_profile_enable = False  
    aoa_graph_enable = True  
    drag_graph_enable = True  
    mass_graph_enable = True  
    thrust_graph_enable = True 

    fast_graphing = False


    # Rocket centric coordinate system definition:
    # Origin point is tip of the nose cone, to physically ground all dimensions.
    # Y axis is the roll axis, with fore being positive and aft being negative.
    # Positive roll is counterclockwise when viewed normal to XZ plane from +Y.
    # X axis is the pitch axis, with starboard being positive and port being negative.
    # Positive pitch is counterclockwise when viewed normal to YZ plane from +X.
    # Z axis is the yaw axis, with dorsal being positive and ventral being negative.
    # Positive yaw is clockwise when viewed normal to XY plane from +Z.

    # constructor
    def __init__(self, name: str, pos: vector, yaw: float, pitch: float, roll: float, v_0: float, ymi: float, pmi: float, rmi: float, cp: vector, cd: float, A: float, cd_s: float, A_s: float, chute_cd: float, chute_A: float,
                 drogue_cd: float, drogue_A: float, cg: vector, dry_mass: float, fuel_mass: float, thrust, t0: float, wind: WindProfile, initDebug: bool, fin: FinSet):
        self.name = name
        # FUNDAMENTAL VECTORS
        self.pos = pos  # m, 3D cartesian position
        self.rot = vec(yaw, pitch, roll)  # rad, orientation
        self.v_0 = v_0  # m/s, initial velocity magnitude
        self.I_0 = vec(ymi, pmi, rmi)  # kg*m^2, Mass moments of inertia
        # COORDINATE SYSTEM
        self.yaw_axis = vec(0,0,1)  # Rocket coordinate system (RocCS) base vectors
        self.pitch_axis = vec(1,0,0)
        self.roll_axis = vec(0,1,0)
        self.global_roll = yp_unit(self.rot.x, self.rot.y)  # RocCS base vectors, transformed to GCS
        self.global_yaw = rotate(yp_unit(self.rot.x, self.rot.y + pi/2), self.rot.z, self.global_roll)
        self.global_pitch = rotate(yp_unit(self.rot.x + pi/2, self.rot.y), self.rot.z, self.global_roll)
        # AERODYNAMIC PROPERTIES
        self.cp = cp  # m, center of pressure position vector (see coordinate system definition above)
        self.cd = cd  # frontal drag coefficient
        self.A = A  # m^2, frontal reference area
        self.cd_s = cd_s  # side drag coefficient (airflow parallel to yaw/pitch axis)
        self.A_s = A_s  # side reference area " "
        self.chute_cd = chute_cd  # parachute drag coefficient
        self.chute_A = chute_A  # m^2, parachute reference area
        self.drogue_cd = drogue_cd  # drogue chute cd
        self.drogue_A = drogue_A  # m^2, drogue chute ref. area
        self.wind = wind  # Wind profile
        self.fin = fin # Fin set
        # MASS PROPERTIES
        self.cg = cg  # m, center of mass position vector (see coordinate system definition above)
        self.dry_mass = dry_mass  # kg
        self.fuel_mass = fuel_mass  # kg
        self.mass = self.dry_mass + self.fuel_mass  # kg, total initial mass
        # PROPULSION PROPERTIES
        self.thrust = thrust  # thrust function of time

        t = 0  # s, thrust time variable
        dt = 0.01  # s, thrust time step
        self.J = 0  # total impulse
        while self.thrust(t) > 0:  # integration of thrust
            self.J += self.thrust(t) * dt
            t += dt
        self.bt = t  # s, motor burn time
        self.t0 = t0  # s, ignition time
        self.t1 = self.t0 + self.bt  # s, burnout time

        # MOMENTUM PROPERTIES
        self.p = self.mass * self.v_0 * yp_unit(self.rot.x, self.rot.y)  # kgm/s, Translational momentum
        self.v = self.p / self.mass  # m/s, Translational velocity
        self.drot = vec(0, 0, 0)  # rad/s, yaw/pitch/roll rate
        self.L = vec(  # Angular momentum, yaw/pitch/roll
            self.I_0.x * self.drot.x,
            self.I_0.y * self.drot.y,
            self.I_0.z * self.drot.z
        )

        # State variables for drogue and main parachute deployment, NOT deployment toggles.
        self.drogue = False
        self.main_chute = False

        # Debug prints
        if initDebug:
            print(f"{self.name} initial momentum: {self.p}kgm/s")
            print(f"{self.name} position: {self.pos}m rotation: {self.rot}rad moments of inertia (YPR): {self.I_0}kgm")
            print(f"{self.name} total impulse: {d2(self.J)}Ns")

        # Graph object declaration

        if FreeRocket.side_profile_enable:
            self.flight_side_graph = graph(title=f"Side Profile of {self.name}", xtitle="Downrange Distance",
                                       ytitle="Altitude (m)", fast=FreeRocket.fast_graphing)
            self.flight_side = gcurve(graph=self.flight_side_graph, color=color.green, label="Side-on Track")

        if FreeRocket.top_profile_enable:
            self.flight_top_graph = graph(title="Top Profile of " + self.name, xtitle="x",
                                          ytitle="z", fast=FreeRocket.fast_graphing)
            self.flight_top = gcurve(graph=self.flight_top_graph, color=color.blue, label="Top-Down Track")

        if FreeRocket.acceleration_graph_enable:
            self.acceleration_graph = graph(title=f"Acceleration of {self.name}", xtitle="t", ytitle="m/s^2", fast=FreeRocket.fast_graphing)
            self.acceleration_x = gcurve(graph=self.acceleration_graph, color=color.red, label="a<sub>x</sub>")  # m/s^2, x-axis acceleration
            self.acceleration_y = gcurve(graph=self.acceleration_graph, color=color.green, label="a<sub>y</sub>")  # m/s^2, y-axis acceleration
            self.acceleration_z = gcurve(graph=self.acceleration_graph, color=color.blue, label="a<sub>z</sub>")  # m/s^2, z-axis acceleration
            self.acceleration_total = gcurve(graph=self.acceleration_graph, color=color.black, label="a")  # m/s^2, total acceleration

        if FreeRocket.position_graph_enable:
            self.position_graph = graph(title=f"Position of {self.name}", xtitle="t", ytitle="m", fast=FreeRocket.fast_graphing)
            self.position_x = gcurve(graph=self.position_graph, color=color.red, label="x")  # m, x position
            self.position_y = gcurve(graph=self.position_graph, color=color.green, label="y")  # m, y position
            self.position_z = gcurve(graph=self.position_graph, color=color.blue, label="z")  # m, z position
            self.position_total = gcurve(graph=self.position_graph, color=color.black, label="Downrange distance")  # m, downrange distance

        if FreeRocket.rotation_graph_enable:
            self.rotation_graph = graph(title=f"Orientation of {self.name}", xtitle="t", ytitle="radians", fast=FreeRocket.fast_graphing)
            self.rotation_yaw = gcurve(graph=self.rotation_graph, color=color.red, label="Yaw")
            self.rotation_pitch = gcurve(graph=self.rotation_graph, color=color.green, label="Pitch")
            self.rotation_roll = gcurve(graph=self.rotation_graph, color=color.blue, label="Roll")

        if FreeRocket.rotation_rate_graph_enable:
            self.rotation_rate_graph = graph(title=f"Rotation rate of {self.name}", xtitle="t", ytitle="radians/s", fast=FreeRocket.fast_graphing)
            self.rotation_rate_yaw = gcurve(graph=self.rotation_rate_graph, color=color.red, label="Yaw rate")
            self.rotation_rate_pitch = gcurve(graph=self.rotation_rate_graph, color=color.green, label="Pitch rate")
            self.rotation_rate_roll = gcurve(graph=self.rotation_rate_graph, color=color.blue, label="Roll rate")
            self.rotation_rate_total = gcurve(graph=self.rotation_rate_graph, color=color.black, label="Total rate")

        if FreeRocket.aoa_graph_enable:
            self.aoa_graph = graph(title=f"AoA of {self.name}", xtitle="t", ytitle="alpha, degrees", fast=FreeRocket.fast_graphing)
            self.aoa = gcurve(graph=self.aoa_graph, color=color.red, label="alpha")  # rad, Angle of attack

        if FreeRocket.moment_graph_enable:
            self.moment_graph = graph(title=f"Moments on {self.name}", xtitle="t", ytitle="Nm", fast=FreeRocket.fast_graphing)
            self.moment_yaw = gcurve(graph=self.moment_graph, color=color.red, label="Yaw moment")
            self.moment_pitch = gcurve(graph=self.moment_graph, color=color.green, label="Pitch moment")
            self.moment_roll = gcurve(graph=self.moment_graph, color=color.blue, label="Roll moment")
            self.moment_total = gcurve(graph=self.moment_graph, color=color.black, label="Total moment")

        if FreeRocket.velocity_graph_enable:
            self.velocity_graph = graph(title=f"Velocity of {self.name}", xtitle="t", ytitle="m/s", fast=FreeRocket.fast_graphing)
            self.velocity_x = gcurve(graph=self.velocity_graph, color=color.red, label="v<sub>x</sub>")
            self.velocity_y = gcurve(graph=self.velocity_graph, color=color.green, label="v<sub>y</sub>")
            self.velocity_z = gcurve(graph=self.velocity_graph, color=color.blue, label="v<sub>z</sub>")
            self.velocity_total = gcurve(graph=self.velocity_graph, color=color.black, label="|v|")

        if FreeRocket.drag_graph_enable:
            self.drag_graph = graph(title=f"Force Due to Drag on {self.name}", xtitle="t", ytitle="N", fast=FreeRocket.fast_graphing)
            self.drag = gcurve(graph=self.drag_graph, color=color.red, label="F<sub>d</sub>")

        if FreeRocket.mass_graph_enable:
            self.mass_graph = graph(title=f"Mass of {self.name}", xtitle="t", ytitle="kg", fast=FreeRocket.fast_graphing)
            self.mass_plot = gcurve(graph=self.mass_graph, color=color.red, label="m")

        if FreeRocket.thrust_graph_enable:
            self.thrust_graph = graph(title=f"Thrust of {self.name}", xtitle="t", ytitle="N", fast=FreeRocket.fast_graphing)
            self.thrust_plot = gcurve(graph=self.thrust_graph, color=color.red, label="Ft")


        self.v_max = 0  # m/s, maximum absolute velocity
        self.v_max_time = 0  # s, max velocity time
        self.y_max = 0  # m, maximum altitude
        self.y_max_time = 0  # s, apogee time
        self.a_max = 0  # m/s^2, maximum acceleration
        self.g_max = 0  # G, maximum acceleration
        self.a_max_time = 0 # s, max accel time
        self.v_ground_hit = 0 # m/s, ground hit velocity
        self.duration = 0 # s, flight duration
        self.q_max = 0  # Pa, maximum dynamic pressure
        self.q_max_time = 0  # s
        self.q_max_speed = 0  # m/s
        self.q_max_accel = 0  # m/s^2
        self.mach_max = 0 # M, maximum Mach number
        self.mach_max_time = 0 # s
        self.mach_max_speed = 0 # m/s
        self.mach_max_altitude = 0 # m
        # end of constructor

    # Drag coefficient estimator
    def cd_alpha(self, alpha):  # Sine interpolation between frontal drag coefficient and side drag coefficient.
        return abs(sin(alpha) * self.cd_s + cos(alpha) * self.cd)

    # Reference area estimator
    def A_alpha(self, alpha):  # Sine interpolation between frontal area and side area.
        return abs(sin(alpha) * self.A_s + cos(alpha) * self.A)

    # Simulation 
    def simulate(self, t, dt):
        # FORCE & MOMENT COMPONENTS

        # BODY DRAG
        heading = yp_unit(self.rot.x, self.rot.y)  # 3d linear unit vector of vehicle orientation
        airflow = self.v + self.wind.wind(self.pos.y)  # 3d linear vector of oncoming airstream (reversed)
        alpha = math.acos(heading.dot(airflow.hat))  # rad, angle of attack
        cd = self.cd_alpha(alpha)  # drag coefficient
        A = self.A_alpha(alpha)  # m^2, reference area
        f_drag = airflow.mag ** 2 * rho(self.pos.y) * cd * A / 2 * -airflow.hat # N, body drag force

        M_drag = f_drag.cross(-self.cg + self.cp)  # Moment generated by drag force

        # FIN LIFT
        M_lift = self.fin.total_lift_vec(self.rot, airflow, self.cg, self.pos.y).cross(self.fin.center - self.cg)

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
        f_grav = vec(0, -g(self.pos.y, self.mass), 0)

        # THRUST VECTOR
        f_thrust = vec(0, 0, 0)
        if self.t0 <= t <= self.t1:
            f_thrust = self.thrust(t - self.t0) * yp_unit(self.rot.x, self.rot.y)

        # TOTAL NET FORCE
        f_net = f_grav + f_thrust + f_drag + f_chute + f_drogue

        # TOTAL NET MOMENT
        M_net = M_drag + M_lift

        self.p = self.p + f_net * dt  # incrementing linear momentum
        self.v = self.p / self.mass  # calculating velocity
        self.pos = self.pos + self.v * dt  # incrementing position

        self.L += M_net * dt  # incrementing angular momentum
        self.drot = vec(self.L.x / self.I_0.x, self.L.y / self.I_0.y, self.L.z / self.I_0.z)  # calculating rotation rate
        self.rot = self.rot + self.drot * dt  # incrementing orientation

        # End of movement during this iteration

        # Remove burned fuel from vehicle mass
        # Assumes fuel mass loss is proportional to thrust
        self.mass = self.mass - f_thrust.mag * dt / self.J * self.fuel_mass

        # Parachute deployment checks

        if self.v.y < 5 and not self.drogue:
            self.drogue = True
            print(f"{self.name} Drogue deployment at T+{d2(t)}s at Altitude:{d2(self.pos.y)}m & Speed:{d2(self.v.mag)}m/s")

        if self.pos.y < 150 and self.drogue and not self.main_chute:
            self.main_chute = True
            print(f"{self.name} Main chute deployment at T+{d2(t)}s at Altitude:{d2(self.pos.y)}m & Speed:{d2(self.v.mag)}m/s")

        mach_n = mag(self.v + self.wind.wind(self.pos.y)) / c(self.pos.y)
        # Graph switched variables
        if FreeRocket.side_profile_enable:
            self.flight_side.plot(sqrt(self.pos.x**2 + self.pos.z**2), self.pos.y)
        if FreeRocket.top_profile_enable:
            self.flight_top.plot(self.pos.x, self.pos.z)
        if FreeRocket.moment_graph_enable:
            self.moment_yaw.plot(t, M_net.x)
            self.moment_pitch.plot(t, M_net.y)
            self.moment_roll.plot(t, M_net.z)
            self.moment_total.plot(t, M_net.mag)
        if FreeRocket.velocity_graph_enable:
            self.velocity_x.plot(t, self.v.x)
            self.velocity_y.plot(t, self.v.y)
            self.velocity_z.plot(t, self.v.z)
            self.velocity_total.plot(t, self.v.mag)
        if FreeRocket.position_graph_enable:
            self.position_x.plot(t, self.pos.x)
            self.position_y.plot(t, self.pos.y)
            self.position_z.plot(t, self.pos.z)
            self.position_total.plot(t, self.pos.mag)
        if FreeRocket.rotation_graph_enable:
            self.rotation_yaw.plot(t, self.rot.x)
            self.rotation_pitch.plot(t, self.rot.y)
            self.rotation_roll.plot(t, self.rot.z)
        if FreeRocket.acceleration_graph_enable:
            self.acceleration_x.plot(t, f_net.x / self.mass)
            self.acceleration_y.plot(t, f_net.y / self.mass)
            self.acceleration_z.plot(t, f_net.z / self.mass)
            self.acceleration_total.plot(t, f_net.mag / self.mass)
        if FreeRocket.rotation_rate_graph_enable:
            self.rotation_rate_yaw.plot(t, self.drot.x)
            self.rotation_rate_pitch.plot(t, self.drot.y)
            self.rotation_rate_roll.plot(t, self.drot.z)
        if FreeRocket.aoa_graph_enable:
            self.aoa.plot(t, alpha/pi*180)
        if FreeRocket.drag_graph_enable:
            self.drag.plot(t, f_drag.mag)
        if FreeRocket.mass_graph_enable:
            self.mass_plot.plot(t, self.mass)
        if FreeRocket.thrust_graph_enable:
            self.thrust_plot.plot(t, f_thrust.mag)

        # Update flight report variables
        self.duration = t
        if self.pos.y > self.y_max:
            self.y_max = self.pos.y
            self.y_max_time = t
        if self.v.mag > self.v_max:
            self.v_max = self.v.mag
            self.v_max_time = t
        a = f_net.mag / self.mass
        if a > self.a_max:
            self.a_max = a
            self.g_max = a / abs(g_0)
            self.a_max_time = t
        if self.pos.y <= 0:
            self.duration = t
            self.v_ground_hit = self.v.mag
        q = 1/2*rho(self.pos.y)*self.v.mag**2
        if self.q_max < q:
            self.q_max = q
            self.q_max_time = t
            self.q_max_speed = self.v.mag
        if self.mach_max < mach_n:
            self.mach_max = mach_n
            self.mach_max_time = t
            self.mach_max_speed = self.v.mag
            self.mach_max_altitude = self.pos.y

        #normalized_airflow = vec(airflow.dot(vec(1,0,0)))
        self.fin.aoa_plot(t, self.rot, airflow, self.cg)

        self.global_roll = yp_unit(self.rot.x, self.rot.y)  # RocCS base vectors, transformed to GCS
        self.global_yaw = rotate(yp_unit(self.rot.x, self.rot.y + pi/2), self.rot.z, self.global_roll)
        self.global_pitch = rotate(yp_unit(self.rot.x + pi/2, self.rot.y), self.rot.z, self.global_roll)

        # end of simulation method

    def flight_report(self):
        print(f"{self.name} Flight Report -------------")
        print(f"Apogee: {d2(self.y_max)}m ({d2(self.y_max*39.37/12)}ft) at T+{d2(self.y_max_time)}s")
        print(f"Maximum speed: {d2(self.v_max)}m/s ({d2(self.v_max*39.37/12/5280*3600)}mph) at T+{d2(self.v_max_time)}s")
        print(f"Maximum acceleration: {d2(self.a_max)}m/s^2 ({d2(self.g_max)}g) at T+{d2(self.a_max_time)}s")
        print(f"Flight Duration: {d2(self.duration)}s")
        print(f"Ground hit velocity: {d2(self.v_ground_hit)}m/s ({d2(self.v_ground_hit*39.37/12)}ft/s), {d2(self.pos.mag)}m downrange")
        print(f"Maximum dynamic pressure: {d2(self.q_max / 1000)}kPa ({d2(self.q_max / 101325 * 14.7)}psig), at {d2(self.q_max_speed)}m/s, at T+{d2(self.q_max_time)}s")
        print(f"Maximum Mach number: {d2(self.mach_max)}M at T+{d2(self.mach_max_time)}s at {d2(self.mach_max_speed)}m/s at {d2(self.mach_max_altitude)}m altitude")

    # end of class def

I435_time_points = [0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4]
I435_thrust_points = [409.4694,447.336,451.1619,450.9657,369.0522,93.0969,76.8123,69.7491,55.2303,36.5913,27.6642,21.4839,17.9523,13.8321,5.4936,2.7468,2.1582]

# thrust interpolation for 38mm 68-16-16 AP-Al motor test fired 3/16/2023
def I435(t):
    if t<=2.4:
        return np.interp(t, I435_time_points, I435_thrust_points)
    else:
        return 0

# thrust function approximation for I240 motor
def I240(t):
    return -77 * (t - 0.2) ** 2 + 300

# thrust function for LR-101 pressure-fed kerolox engine
def LR101(t):
    if t < 11:
        return 830 * 4.448
    else:
        return 0

# alpha phoenix fin coefficient of lift function
def cl(aoa: float):
    return 2*pi*aoa # using NASA approx for thin subsonic airfoils

# Beginning of actual program execution

AlphaPhoenixFins = dict(num_fins=4, center=vec(0,-0.152,0), pos=vec(0,-0.162,0), planform=0.005, stall_angle=10*pi/180, ac_span=0.05, cl_pass=cl)
fin_1 = FinSet(**AlphaPhoenixFins)

FAR_wind = dict(name="FAR", mu=3, sigma=2, angle_mu=pi/4, angle_sigma=pi/8, step=100, print_debug=False)
wind_1 = WindProfile(**FAR_wind)

# Alpha Phoenix on I-300
AlphaPhoenix = dict(name="Alpha Phoenix", pos=vec(0,1,0), yaw=0, pitch=90*pi/180, roll=0, v_0=5, ymi=0.0715, pmi=0.0715, rmi=0.0012, cp=vec(0,-0.152,0), cd=0.6, A=0.015, cd_s=1.5, A_s = 0.05, chute_cd=0.8, chute_A=0.4, drogue_cd=0.8, drogue_A=0.1, cg=vec(0,-0.142,0), dry_mass=1.1, fuel_mass=0.220, thrust=I240, t0=0, wind=wind_1, initDebug=True, fin=fin_1)

# Theseus on LR-101
TheseusFins= dict(num_fins=4, center=vec(0,-5,0), pos=vec(0,-5.2,0), planform=0.0258, stall_angle=10*pi/180, ac_span=0.165, cl_pass=cl)
fin_theseus = FinSet(**TheseusFins)

Theseus = dict(name="Theseus", pos=vec(0,1,0), yaw=0, pitch=90*pi/180, roll=0, v_0=5, ymi=0.0715*100, pmi=0.0715*100, rmi=0.0012*100, cp=vec(0,-4.5,0), cd=0.425, A=0.0411, cd_s=1, A_s=0.5, chute_cd=0.8, chute_A=1, drogue_cd=0.8, drogue_A=0.1, cg=vec(0,-4.5+8/39.37,0), dry_mass=200/2.204, fuel_mass=46.345/2.204, thrust=LR101, t0=0, wind=wind_1, initDebug=True, fin=fin_theseus)

payload = FreeRocket(**Theseus)

graphics_3D = False
if graphics_3D:
    #payload_track = sphere(pos=payload.pos, size=vec(1,1,1), color=color.red, make_trail=True)
    yaw_arrow = arrow(pos=vec(0,0,0), axis=payload.global_yaw, color=color.red)
    pitch_arrow = arrow(pos=vec(0,0,0), axis=payload.global_pitch, color=color.green)
    roll_arrow = arrow(pos=vec(0,0,0), axis=payload.global_pitch, color=color.blue)
    #v_arrow = arrow(pos=vec(0,0,0), axis=payload.v * 10, color=color.black)

time = 0
dtime = 1 / 20

print(f"Time step: {dtime}")
print("----BEGIN SIMULATION----")
while time < 120 and not payload.drogue:
    #sleep(1/8)
    payload.simulate(time, dtime)
    if graphics_3D:
        #payload_track.pos = payload.pos
        #yaw_arrow.pos = payload.pos
        yaw_arrow.axis = payload.global_yaw
        #pitch_arrow.pos = payload.pos
        pitch_arrow.axis = payload.global_pitch
        #roll_arrow.pos = payload.pos
        roll_arrow.axis = payload.global_roll

    time += dtime
print("-----END SIMULATION-----")

payload.flight_report()

