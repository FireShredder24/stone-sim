# By John Nguyen <knightofthealtar64@gmail.com>
# Python 3.10.11
# MIT License

# Regulated pressure-fed fluids system COPV pressure calculator

# This program aims to simulate the flow of pressurant gas through the dome regulators
# of a pressure-fed rocket fluids system in order to determine the initial amount of 
# pressure required to maintain pressure throughout the entire burn.

import numpy as np
from vpython import *

# Universal constants

R_btu = 1.985875279009 # BTU / (lbmol * R)
R_inlb = R_btu * 9338.0318994166 # in-lbf / (lbmol * R)

# Clausius-Clapeyron equation for vapor pressure
def vapor_pressure(h_vap, M, boil, T):
    return 14.7 * np.exp(-h_vap * M / R_btu * (1 / T - 1 / boil))

# System definitions

copv_vol = 416 # cu. in.
ox_tank_vol = 980 # cu. in.
fuel_tank_vol = 814 # cu. in.

ox_dr_press = 600 + 14.7 # psia, ox side dome regulator pressure setting
fuel_dr_press = 565 + 14.7 # psia, fuel side dome regulator pressure setting

T_amb = 530 # R, ambient temperature

copv_press = 3000
copv_T = T_amb # R, initial COPV temperature

ox_ul_press = 600 + 14.7 # psia, initial ox tank ullage gas pressure
fuel_ul_press = 565 + 14.7 # psia, initial fuel tank ullage gas pressure
ox_total_press = ox_ul_press # psia, initial total pressure
fuel_total_press = fuel_ul_press # psia, initial total pressure

ox_press_drop = 58 # psia, pressure drop over piping between dome regulator and tank
fuel_press_drop = 30 # psia, pressure drop over piping between dome regulator and tank

ox_ul_vol = 100 # cu. in., initial gas volume in ox tank
fuel_ul_vol = 100 # cu. in., initial gas volume in fuel tank

ox_vol_rate = 0.0376 * 12**3 # in^3/s, liquid oxidizer volume flow rate
fuel_vol_rate = 0.0291 * 12**3 # in^3/s, liquid fuel volume flow rate

# Liquid oxidizer properties

ox_name = "Liquid Oxygen"

ox_M = 32 # lbm / lbmol, molar mass of oxidizer
ox_T_boil = 162.4 # R, oxidizer boiling point
ox_T = 200 # R, oxidizer initial temperature
ox_rho_l = 71.16 / 12**3 # lb / in^3
ox_heat_vap = 91.7 # BTU / lbm
ox_vap_press = vapor_pressure(ox_heat_vap, ox_M, ox_T_boil, ox_T) # psia, ox vapor partial pressure
ox_total_press += ox_vap_press # psia, ox tank total pressure
ox_c_l = 0.219 # BTU / lbm / R
ox_cp = 0.22 # BTU / lbm / R, specific heat at constant pressure
ox_cv = 0.16 # BTU / lbm / R, specific heat at constant volume

# Liquid fuel properties

fuel_name = "Kerosene"

fuel_M = 177 # lbm / lbmol, molar mass of fuel
fuel_T = T_amb # R, fuel initial temperature
fuel_T_boil = 770 # R, fuel boiling point
fuel_rho_l = 51.1 / 12**3 # lb / in^3, liquid phase density
fuel_heat_vap = 108 # BTU / lbm, heat of vaporization
fuel_vap_press = vapor_pressure(fuel_heat_vap, fuel_M, fuel_T_boil, fuel_T) # psia, fuel vapor partial pressure
fuel_total_press += fuel_vap_press # psia, fuel tank total pressure
fuel_c_l = 0.48 # BTU / lbm / R
# values for kerosene unknown.  Using same as nitrogen since kerosene unlikely to vaporize sigificantly
fuel_cp = 0.25 # BTU / lbm / R, specific heat at constant pressure
fuel_cv = 0.18 # BTU / lbm / R, specific heat at constant volume
# Ullage gas properties

ul_name = "Nitrogen"

ul_M = 28 # lbm / lbmol, molar mass of pressurant gas
ul_T_boil = 138.7 # R, ullage gas boiling point
ul_T = T_amb # R, ullage gas starting temperature
ul_rho_l = 50 / 12**3 # lbm / in^3, liquid phase density
ul_heat_vap = 86 # BTU / lbm, heat of vaporization
ul_cp = 0.25 # BTU / lbm / R, specific heat at constant pressure
ul_cv = 0.18 # BTU / lbm / R, specific heat at constant volume

ox_vapor_mass = ox_M * ox_vap_press * ox_ul_vol / R_inlb / ox_T
ox_ul_mass = ul_M * ox_ul_press * ox_ul_vol / R_inlb / ox_T

fuel_vapor_mass = fuel_M * fuel_vap_press * fuel_ul_vol / R_inlb / fuel_T
fuel_ul_mass = ul_M * fuel_ul_press * fuel_ul_vol / R_inlb / fuel_T

copv_ul_mass = ul_M * copv_press * copv_vol / R_inlb / copv_T

print(f"Regulated pressure-fed rocket feed system simulator")
print(f"COPV: {copv_vol} in^3 of {ul_name} at {copv_press} psia and {ul_T} R")
print(f"Oxidizer: {ox_ul_vol}/{ox_tank_vol} in^3 of {ox_name} at {ox_T} R with initial pressure {ox_total_press:.2f} psia")
print(f"Fuel: {fuel_ul_vol}/{fuel_tank_vol} in^3 of {fuel_name} at {fuel_T} R with initial pressure {fuel_total_press:.2f} psia")
print(f"Oxidizer flow: {ox_vol_rate:.2f} in^3/s at {ox_dr_press} psia")
print(f"Fuel flow: {fuel_vol_rate:.2f} in^3/s at {fuel_dr_press} psia")

press_graph = graph(title="Pressure", fast=False, xtitle="s", ytitle="psia")
copv_press_curve = gcurve(graph=press_graph, label="COPV", color=color.green)
ox_press_curve = gcurve(graph=press_graph, label="OX", color=color.blue)
ox_vap_press_curve = gcurve(graph=press_graph, label="OX VAP", color=color.purple)
fuel_press_curve = gcurve(graph=press_graph, label="FUEL", color=color.red)
fuel_vap_press_curve = gcurve(graph=press_graph, label="FUEL VAP", color=color.orange)

temp_graph = graph(title="Temperature", fast=False, xtitle="s", ytitle="R")
copv_T_curve = gcurve(graph=temp_graph, label="COPV", color=color.green)
ox_T_curve = gcurve(graph=temp_graph, label="OX", color=color.blue)
ox_dr_T_curve = gcurve(graph=temp_graph, label="OX DR", color=color.purple)
fuel_T_curve = gcurve(graph=temp_graph, label="FUEL", color=color.red)
fuel_dr_T_curve = gcurve(graph=temp_graph, label="FUEL DR", color=color.orange)

volume_graph = graph(title="Volume", fast=False, xtitle="s", ytitle="in^3")
ox_vol_curve = gcurve(graph=volume_graph, label="OX", color=color.blue)
fuel_vol_curve = gcurve(graph=volume_graph, label="FUEL", color=color.red)

mass_graph = graph(title="Mass", fast=False, xtitle="s", ytitle="lbm")
copv_mass_curve = gcurve(graph=mass_graph, label="COPV", color=color.green)
ox_ul_mass_curve = gcurve(graph=mass_graph, label="OX UL", color=color.blue)
ox_vap_mass_curve = gcurve(graph=mass_graph, label="OX VAP", color=color.purple)
fuel_ul_mass_curve = gcurve(graph=mass_graph, label="FUEL UL", color=color.red)
fuel_vap_mass_curve = gcurve(graph=mass_graph, label="FUEL VAP", color=color.orange)

t = 0
dt = 0.025

print("Begin simulation")
while ox_ul_vol <= ox_tank_vol and fuel_ul_vol <= fuel_tank_vol and t < 30:
    # Calculating vapor pressures of liquid propellants.
    # We assume that the relatively violent pressurization will cause instant saturation.
    # Note that the temperature here corresponds to the temperature of the gas above liquid.
    # This is added to the total pressure of the tank
    ox_vap_press = vapor_pressure(ox_heat_vap, ox_M, ox_T_boil, ox_T)
    fuel_vap_press = vapor_pressure(fuel_heat_vap, fuel_M, fuel_T_boil, fuel_T)

    # First we need to calculate the temperature at the dome regulator outlets
    ox_dr_T = copv_T * (ox_dr_press / copv_press) ** (1-ul_cv/ul_cp)
    fuel_dr_T = copv_T * (fuel_dr_press / copv_press) ** (1-ul_cv/ul_cp)

    # Then using this we determine the mass of ullage gas being added to the tanks
    if ox_dr_press > ox_total_press:
        ox_new_mass = ul_M * min(copv_press,ox_dr_press) * ox_vol_rate / R_inlb / ox_dr_T * dt
    else:
        ox_new_mass = 0
    if fuel_dr_press > fuel_total_press:
        fuel_new_mass = ul_M * min(copv_press, fuel_dr_press) * fuel_vol_rate / R_inlb / fuel_dr_T * dt
    else:
        fuel_new_mass = 0

    # Then, due to pressure drop over the pipework, the temperature of the gas entering the tanks will be lower
    ox_dr_T = ox_dr_T * ((ox_dr_press - ox_press_drop) / ox_dr_press) ** (1-ul_cv/ul_cp)
    fuel_dr_T = fuel_dr_T * ((fuel_dr_press - fuel_press_drop) / fuel_dr_press) ** (1-ul_cv/ul_cp)

    # Then we take the existing mass within the tanks and determine the final temperature

    ox_T = (ox_vapor_mass*ox_cv*ox_T + ox_ul_mass*ul_cv*ox_T + ox_new_mass*ul_cv*ox_dr_T) / (ox_vapor_mass*ox_cv + (ox_ul_mass + ox_new_mass)*ul_cv)

    fuel_T = (fuel_vapor_mass*fuel_cv*fuel_T + fuel_ul_mass*ul_cv*fuel_T + fuel_new_mass*ul_cv*fuel_dr_T) / (fuel_vapor_mass*fuel_cv + (fuel_ul_mass + fuel_new_mass)*ul_cv)

    # Then we find final pressure
    ox_ul_press = (ox_ul_mass + ox_new_mass) / ul_M * R_inlb * ox_T / ox_ul_vol
    ox_total_press = ox_ul_press + ox_vap_press

    fuel_ul_press = (fuel_ul_mass + fuel_new_mass) / ul_M * R_inlb * fuel_T / fuel_ul_vol
    fuel_total_press = fuel_ul_press + fuel_vap_press

    # And increase the ullage volume available
    ox_ul_vol += ox_vol_rate * dt
    fuel_ul_vol += fuel_vol_rate * dt

    # And now update the COPV pressure and temperature
    copv_ul_mass -= ox_new_mass + fuel_new_mass
    copv_press_old = copv_press
    copv_press = copv_ul_mass / ul_M * R_inlb * copv_T / copv_vol
    copv_T = copv_T*(copv_press / copv_press_old) ** (1-ul_cv/ul_cp)

    # Update vapor and ullage mass

    ox_vapor_mass = ox_M * ox_vap_press * ox_ul_vol / R_inlb / ox_T
    ox_ul_mass += ox_new_mass

    fuel_vapor_mass = fuel_M * fuel_vap_press * fuel_ul_vol / R_inlb / fuel_T
    fuel_ul_mass += fuel_new_mass

    # Printing for posterity
    print(f"t={t:.2f}, copv: {copv_press:.0f} psia & {copv_T:.0f} R, ox at {ox_total_press:.0f} psia & {ox_T:.0f} R, fuel at {fuel_total_press:.0f} psia & {fuel_T:.0f} R")

    # Graphing
    copv_press_curve.plot(t, copv_press)
    copv_T_curve.plot(t, copv_T)
    copv_mass_curve.plot(t, copv_ul_mass)
    ox_press_curve.plot(t, ox_total_press)
    ox_vap_press_curve.plot(t, ox_vap_press)
    ox_T_curve.plot(t, ox_T)
    ox_dr_T_curve.plot(t, ox_dr_T)
    ox_ul_mass_curve.plot(t, ox_ul_mass + ox_new_mass)
    ox_vap_mass_curve.plot(t, ox_vapor_mass)
    fuel_press_curve.plot(t, fuel_total_press)
    fuel_vap_press_curve.plot(t, fuel_vap_press)
    fuel_T_curve.plot(t, fuel_T)
    fuel_dr_T_curve.plot(t, fuel_dr_T)
    fuel_ul_mass_curve.plot(t, fuel_ul_mass + fuel_new_mass)
    fuel_vap_mass_curve.plot(t, fuel_vapor_mass)
    ox_vol_curve.plot(t, ox_ul_vol)
    fuel_vol_curve.plot(t, fuel_ul_vol)



    # and of course increment time
    t += dt
