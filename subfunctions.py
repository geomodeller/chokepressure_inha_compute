## 1.find pressure kick of bottom
## kick은 annulus에 위치해 있음 따라서 애눌러스의 부피 부분을 구하고 킥의 부피를 가정하면 kick의 depth를 구한다
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from math import log, nan

def inside_capacity(openhole, pipe_OD):
    ID_c = (openhole**2- pipe_OD**2)/ 1029.4 
    return ID_c

def find_depth_between(ID_c, volume):
    f_depth_between = volume / ID_c
    return f_depth_between

def Pressure_kick_bottom(kick_height_from_bh, bhp, m_rho):
    P_k_b = bhp - 0.052 * kick_height_from_bh * m_rho 
    return P_k_b

def compute_gas_density(z, R, T, M, pressure):
    return pressure*101.3/14.7/(z*R*T)*M/1000*8.33

def pressure_kick_z(volume, z, n, R, T):
    p_k_z = 0.145 * (z * n * R * T )/(0.1589 *volume)
    return p_k_z

def volume_kick_z(pressure, z, n, R, T):
    return 0.145 * (z * n * R * T )/(0.1589 *pressure)
    

def P_hy_between(k_rho, f_depth_between):
    P_hy_bw = 0.052 * k_rho * f_depth_between / 2 
    return P_hy_bw

def P_kick_middle_r(P_k_b, P_hy_bw):
    P_k_m_r = P_k_b - P_hy_bw
    return P_k_m_r 

def find_matching_volume(P_k_b, P_hy_bw, P_k_m_r, p_k_z, k_rho, depth, m_rho,T, volume_initial, z, ID_c, n, R, f_depth_between):
    step = 0.00001
    volume = volume_initial
    
    while abs(P_k_m_r - p_k_z) > 0.01:
        volume += step
        f_depth_between = find_depth_between(ID_c, volume)
        P_hy_bw = P_hy_between(k_rho, f_depth_between)
        P_k_b = Pressure_kick_bottom(depth, f_depth_between, m_rho, k_rho)
        P_k_m_r = P_kick_middle_r(P_k_b, P_hy_bw)
        p_k_z = pressure_kick_z(volume, z, n, R, T)
    return volume

## 4. find pressure kick of top
def pressure_kick_top(cal_P_k_m_r, cal_P_hy_bw):
    p_k_t = cal_P_k_m_r - cal_P_hy_bw
    return p_k_t

def velocity(q, openhole, pipe_OD):
    vel = q / 2.448 / (openhole**2 - pipe_OD**2)
    return vel

def apparent_viscosity(vel, openhole, pipe_OD, tau_y, mu_p):
    mu_a = mu_p + 5*tau_y*(openhole - pipe_OD)/vel
    return mu_a
 
def reynolds(m_rho, vel, openhole, pipe_OD, mu_a):
    Re = 928* 0.816 * m_rho*vel*(openhole - pipe_OD)/mu_a
    return Re

def friction_factor(Re, rel_roughness):
    if Re>4000:
        # when it's turbulent
        A_0 = -0.79638*log(rel_roughness/8.208 + 7.3357/Re)
        A_1 = Re*rel_roughness + 9.3120665*A_0
        f = ((8.128943 + A_1)/(8.128943*A_0 + 0.86859209*A_1*log(A_1/3.7099535/Re)))**2
    elif Re < 2100:
        # When it's laminar 
        f = 64/Re
    else:
        f = nan
    return f/4    

def pressure_gradient_laminar(mu_p, tau_y, vel, openhole, pipe_OD):
    return mu_p*vel/(1000* (openhole- pipe_OD)**2) + tau_y/(200 *(openhole - pipe_OD)) 

def pressure_gradient_turbulent(f,m_rho, vel, openhole, pipe_OD):
    return f*m_rho*vel**2/21.1/(openhole - pipe_OD)

def Pressure_friction(dP_per_dL, depth):
    return dP_per_dL* (10000 - depth)

## 1. shut in
def choke_pressure1(bhp, m_rho, cal_f_depth_between, k_rho):
    ch_pressure1 = bhp - 0.052 * m_rho * (10000 - cal_f_depth_between) - 0.052 * k_rho * cal_f_depth_between
    return ch_pressure1
   
def cal_find_depth_between(ID_c, v_k_x):
    cal_f_depth_between = v_k_x / ID_c 
    return cal_f_depth_between

def cal_P_hy_between(k_rho, cal_f_depth_between):
    cal_P_hy_bw = 0.052 * k_rho * cal_f_depth_between 
    return cal_P_hy_bw

def cal_Pressure_kick_bottom(depth, cal_f_depth_between, m_rho, k_rho, sicp):
    cal_P_k_b = 0.052 * (depth - cal_f_depth_between) * m_rho + 0.052 * cal_f_depth_between * k_rho + sicp
    return cal_P_k_b

def cal_P_kick_middle_r(cal_P_k_b, cal_P_hy_bw):
    cal_P_k_m_r = cal_P_k_b - (cal_P_hy_bw)/ 2
    return cal_P_k_m_r 
    
def pressure_kick_top(cal_P_k_m_r, cal_P_hy_bw):
    p_k_t = cal_P_k_m_r - cal_P_hy_bw/2
    return p_k_t
    
def volume_kick(new_volume, bhp, p_k_t):
    v_k_x = new_volume * bhp / p_k_t
    return v_k_x

def choke_pressure2(p_k_t_new, m_rho, depth, dP_per_dL, r_cal_f_depth_between_new):
    ch_pressure2 = p_k_t_new - 0.052 * m_rho * (depth - r_cal_f_depth_between_new) - dP_per_dL * (10000-depth)
    return ch_pressure2

    
def choke_pressure3(p_k_t_new, m_rho, depth, dP_per_dL, r_cal_f_depth_between_new):
    ch_pressure3 = p_k_t_new - 0.052 * m_rho * (depth - r_cal_f_depth_between_new) - dP_per_dL * (10000-depth)
    return ch_pressure3