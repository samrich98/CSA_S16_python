from handcalcs.decorator import handcalc
from math import sqrt, pi, cos, asin
from IPython.display import display, Math # to manually rendering equations in latex when handcalcs is insufficient
import forallpeople # units library - from https://github.com/connorferster/forallpeople
forallpeople.environment('structural', top_level=True) # get structural units

# define forallpeople units
MPa = MPa  # type: ignore
inch = inch  # type: ignore
mm = mm  # type: ignore
N = N  # type: ignore

#=== 13.1 Resistance factors ===#
phi = 0.9 # structural steel (yield strength)
phi_u = 0.75 # structural steel (ultimate strength)
phi_b = 0.80 # bolts
phi_br = 0.80 # bearing of bolts on steel
phi_bi = 0.80 # beam web bearing, interior
phi_be = 0.75 # beam web bearing, end

#=== 13.2 Axial tension ===#
@handcalc(jupyter_display = False, precision = 2)
def T_r_y_func(A_g, F_y):
    """Calculates the tensile resistance per clause 13.2. a).
    
    Parameters:
    - A_g = Gross cross-section area
    - F_y = Member yield strength"""
    
    T_r = phi * A_g * F_y
    return T_r

@handcalc(jupyter_display = False, precision = 2)
def T_r_u_func(A_ne, F_u):
    """Calculates the tensile resistance per clause 13.2. c)
    
    Parameters:
    - A_ne = Effective net area from clause 12.3.3
    - F_u = Member ultimate strength"""
    
    T_r = phi_u * A_ne * F_u
    return T_r

#=== 13.3 Axial compression ===#
@handcalc(jupyter_display = False, precision = 2)
def C_r_func(A, F_y, F_e, n = 1.34):
    """The factored axial compressive resistance per clause 13.3.1.1
    
    Parameters:
    - A = Cross-section area
    - F_y = Member yield strength
    - F_e = Effective yield strength
    - n = 1.34 by default, can be changed to 2.24 as permitted in clause 13.3.1.1"""
    
    lamb = sqrt(F_y / F_e)
    
    C_r = (phi * A * F_y) / (1 + lamb ** (2 * n))**(1 / n)
    
    return C_r

@handcalc(jupyter_display = False, precision = 2)
def F_e_func(K, L, r, E = 200000 * MPa):
    """Gives the effective strength of a member under a compressive axial load per clause 13.3.1.2
    
    Parameters:
    - K = Effective length factor
    - L = Unbraced length
    - r = Radius of gyration
    - E = Modulus of elasticity (200,000 MPa by default)"""
    
    F_e = (pi ** 2 * E) / ((K * L) / r) ** 2
    
    return F_e

#=== 13.4 Shear ===#
@handcalc(jupyter_display = False, precision = 2)
def V_r_unstiffened(A_w, h, w, F_y, MPa = MPa):
    """Calculates the shear resistance of a web of a flexureal member with two flanges per clause 13.4
    
    Parameters:
    - A_w = Shear area (dw for rolled shapes, hw for girders, and 2ht for rectangular HSS)
    - h = Web height
    - w = Web width
    - F_y = Member yield strength"""
    
    if h / w <= 1014 / sqrt(F_y): F_s = 0.66 * F_y
    elif h / w <= 1435 / sqrt(F_y): Fy = F_y / MPa; F_s = ((670 * sqrt(Fy)) / (h / w)) * MPa
    elif h / w > 1435 / sqrt(F_y): Fy = F_y / MPa; F_s = (961200 / (h / w) ** 2) * MPa
    
    V_r = phi * A_w * F_s
    
    return V_r

#=== 13.5 Bending ===#
@handcalc(jupyter_display = False, precision = 2)
def M_r_func(Z, F_y):
    """Calculates the bending resistance per clause 13.5.
    Parameters:
    - Z = Section modulus
    - F_y = Member yield strength
    - L = Span length"""
    
    M_r = phi * Z * F_y

    return M_r

#=== 13.13 Block Shear ===#
@handcalc(jupyter_display = False, precision = 2)
# long
def block_shear(U_t, A_n, A_gv, F_y, F_u, MPa = MPa):
    """Calculates the block shear resistance of a tension member per clause 13.11.
    
    Parameters:
    - U_t = Efficiecy factor and U_t = 1.0 for symmetrical blocks or failure patters and concentric loading, 
    otherwise refer to clause 13.11 a)
    - A_n = net area in tension as specified in clause 12
    - A_gv = gross area in shear
    - F_y = Member yield strength
    - F_u = Member ultimate strength"""
    
    if F_y > 460 * MPa: T_r = phi_u * (U_t * A_n * F_u + 0.6 * A_gv * F_y)
    elif F_y <= 460 * MPa: T_r = phi_u * (U_t * A_n * F_u + 0.6 * A_gv * (F_y + F_u) / 2)
        
    return T_r

#=== 21 Connections ===#
## 21.3 Moment-connected members ##
def stiffener_check_comp_flange(F_yc, w_c, t_b, h_c = 0, End = False, Slender = False):
    """Gives the bearing capacity of the column opposite to the beam compression flange in a moment-connnected 
    member per clause 21.3.
    
    Parameters:
    - F_yc = Column yield strength
    - w_c = Thickness of the column web
    - t_b = Thickness of the beam flange
    - h_c = Height of the column web
    - End = Boolean value indicating if the beam is near the top of the column
    - Slender = Boolean value. True if the column web is class 3 or 4"""
    
    if End == False and Slender == False:
        B_r = stiffener_check_comp_flange_compact(F_yc, w_c, t_b)
    elif End == True and Slender == False:
        B_r = stiffener_check_comp_flange_compact_end(F_yc, w_c, t_b)
    elif End == False and Slender == True:
        B_r = stiffener_check_comp_flange_slender(F_yc, w_c, t_b, h_c)
    else:
        B_r = stiffener_check_comp_flange_slender_end(F_yc, w_c, t_b, h_c)
    
    return B_r

@handcalc(jupyter_display = False, precision = 2)
def stiffener_check_comp_flange_compact(F_yc, w_c, t_b):
    """Class 1 or 2 column web with a beam far from the column end."""
    
    B_r = phi_bi * w_c * (t_b + 10)*F_yc
    return B_r

@handcalc(jupyter_display = False, precision = 2)
def stiffener_check_comp_flange_compact_end(F_yc, w_c, t_b):
    """Class 1 or 2 column web with a beam near to the column end."""
        
    phi_bi = 0.75
    B_r = phi_bi * w_c * (t_b + 4)*F_yc
    return B_r

@handcalc(jupyter_display = False, precision = 2)
def stiffener_check_comp_flange_slender(t_c, w_c, t_b, d_c, N = N, mm = mm):
    """Class 1 or 2 column web with a beam far from the column end."""
    
    B_r = (640000 * phi_bi * w_c * (t_b + 10 * t_c))/((d_c - 2 * t_c)/w_c)**2 * N / mm ** 2
    return B_r

@handcalc(jupyter_display = False, precision = 2)
def stiffener_check_comp_flange_slender_end(t_c, w_c, t_b, d_c, N = N, mm = mm):
    """Class 3 or 4 column web with a beam near to the column end."""
    
    phi_bi = 0.75
    B_r = (640000 * phi_bi * w_c * (t_b + 4 * t_c))/((d_c - 2 * t_c)/w_c)**2 * N / mm ** 2
    return B_r

@handcalc(jupyter_display = False, precision = 2)
def stiffener_check_tens_flange(t_c, F_yc):
    """Gives the bearing capacity of the column opposite to the beam tension flange in a moment-connnected 
    member per clause 21.3.
    
    Parameters:
    - F_yc = Column yield strength
    - t_c = Thickness of the column flange"""
    
    T_r = 7 * phi * t_c ** 2 * F_yc
    
    return T_r

def stiffeners_check(B_r, T_r, M_f, d_b, b_c, t_c):
    """Checks if stiffeners are required per clause 21.3 and returns required force developed by the 
    stiffeners.
    
    Parameters:
    - B_r = Bearing resistance of the column flange
    - T_r = Tensile resistance of the column flange
    - M_f = Factored moment acting on the column flange
    - d_b = Depth of the beam
    - b_c = Width of the beam
    - t_c = thickness of the column"""

    if B_r < M_f / d_b: result_1 = "ok" # "ok" if the capacity is greater than the demand
    else: result_1 = "not ok" # otherwise not ok
        
    display(Math(fr"B_r={B_r:.1f} < \frac{{M_f}}{{d_b}}=\frac{{{M_f}}}{{{d_b}}}={M_f / d_b:.1f}" + 
                 fr"\quad \textbf{{{result_1}}}")) # render math

    if T_r < M_f / d_b: result_2 = "ok" # "ok" if the capacity is greater than the demand
    else: result_2 = "not ok" # otherwise not ok
        
    display(Math(fr"T_r={T_r:.1f} < \frac{{M_f}}{{d_b}}=\frac{{{M_f}}}{{{d_b}}}={M_f / d_b:.1f}" + 
                 fr"\quad \textbf{{{result_2}}}")) # render math
    
    if result_1 == "ok" and result_2 == "ok":
        display(Math(fr"\text{{Web stiffeners are not required}}"))
        return 0
    else:
        display(Math(fr"\text{{Web stiffeners are required. Stiffeners must develop a force }}F_{{st}}  " + 
                     fr"\text{{ equal to the maximum of:}}"))
        F_st = F_st_calc(B_r, T_r, M_f, d_b, b_c, t_c) # return the maximum force on the stiffner plates
        return F_st

@handcalc(jupyter_display = False, precision = 2)
def F_st_calc(B_r, T_r, M_f, d_b, b_c, t_c):
    """This function calculated and renders the maximum force that must be developed by the stiffener plates in
    clause 21.3"""
    
    F_st_1 = (M_f / d_b) - B_r
    F_st_2 = (M_f / d_b) - T_r
    F_st_3 = 0.25 * (M_f / d_b)
    F_st_4 = (b_c - 4 * t_c) / b_c * (M_f / d_b)
    
    return max(F_st_1, F_st_2, F_st_3, F_st_4)

#=== 22 Design and Detailing of bolted connections ===#
def bolt_spacing(d_bt,
                 s_edge,
                 s_1, 
                 s, 
                 t_p, 
                 bolt_lines = 3, 
                 edge_condition = "rolled", 
                 end_condition = "sheared", 
                 mm = mm):
    """Checks the minimum picth, minimum and maximum edge distance, and minimum end distance in accordance with
    clause 22.3.
    
    Parameters:
    - d_bt = Diamter of the bold
    - s_edge = Edge distance
    - s_1 = End distance
    - s = Bolt spacing
    - t_p = Thickness of the outside connected element
    - bolt_lines = Number of bolts parallel to the line of action of the force (assumed as greater than 2)
    - edge_condition = Either sheared or rolled, depending on the member edge condition (assumed as rolled)
    - end_condition = Either sheared or rolled, depending on the member edge condition (assumed as sheared)"""
    
    # Minimum pitch clause 22.3.1
    s_min = 2.7 * d_bt
    if s >= s_min: result = "ok"
    else: result = "not ok"
    
    display(Math(r"\hspace{0pt}\text{Minimum pitch:}"))
    display(Math(fr" s = {s:.1f} \geq 2.7d_{{bt}} = 2.7·{d_bt:.1f}={(2.7 * d_bt):.1f} " + 
                 fr"\quad \textbf{{{result}}}\tag*{{cl 22.3.1}}"))
    
    #Minimum edge distance 22.3.2
    edge_min = table_5(d_bt, edge_condition)
    if s_edge >= edge_min: result = "ok"
    else: result = "not ok"

    display(Math(fr"\text{{Minimum edge distance from CSA S16 Table 5 is }}{edge_min:.1f}:"))
    display(Math(fr" s_{{edge}} = {s_edge:.1f} \geq {edge_min:.1f} \quad \textbf{{{result}}}\tag*{{cl 22.3.2}}"))

    #Minimum edge distance 22.3.3
    edge_max = min(12 * t_p, 150 * mm)
    if s_edge <= edge_max: result = "ok"
    else: result = "not ok"
    
    display(Math(fr"\text{{Maximum edge distance: ({edge_condition} edge)}}"))
    display(Math(fr" s_{{edge}} = {s_edge:.1f} \leq min(12t_p, 150 mm) = min(12{t_p:.1f}, 150 mm) = " + 
                 fr"{edge_max:.1f} \quad \textbf{{{result}}}\tag*{{cl 22.3.3}}"))
        
    # Minimum end distance 22.3.4
    if bolt_lines <= 2:
        s_1_min = 1.5 * d_bt
        if s_1 >= s_1_min: result = "ok"
        else: result = "not ok"
            
        display(Math(fr"\text{{Minimum end distance (≤ 2 bolts parallel to the direction of " + 
                     fr"the load): \text{{{end_condition} edge}}"))
        display(Math(fr" s = {s_1:.1f} \geq 1.5d_{{bt}} = 1.5·{d_bt:.1f}={(1.5 * d_bt):.1f} " + 
                     fr"\quad \textbf{{{result}}}\tag*{{cl 22.3.4}}"))
    else:
        s_1_min = table_5(d_bt, end_condition)
        if s_1 >= s_1_min: result = "ok"
        else: result = "not ok"
        
        display(Math(fr"\text{{Minimum end distance (> 2 bolts parallel to the direction of " + 
                     fr"the load) from CSA S16 Table 5 is }}{s_1_min:.1f}: \text{{{end_condition} edge}}"))
        display(Math(fr" s_{{1}} = {s_1:.1f} \geq {s_1_min:.1f} \quad \textbf{{{result}}}\tag*{{cl 22.3.4}}"))
        
        
#=== 27 Seismic Design ===#
## Probable Yield Stress ##
def R_y_func(F_y, HSS = False):
    """Gives the probable yield stress factor, R_y per clause 27.1.7. The minimum probable yield stress is 
    385 MPa for all sections, except HSS.  HSS (HSS = True) have a probable yield stress of 460 MPa or greater.
    
    Parameters:
    - F_y = Member yield strength
    - HSS = Boolean value, if True than this function changes the minimum probable yield stress to 460 MPa."""
    
    if HSS:
        R_y = R_y_calc(F_y, 460 * MPa)
    else:
        R_y = R_y_calc(F_y, 385 * MPa)
        
    return R_y

@handcalc(jupyter_display = False, precision = 2)
def R_y_calc(F_y, F_y_min):
    """Gives the probable yield stress factor, R_y per clause 27.1.7 based on the yield strength and minimum
    allowable probable yield strength.
    
    Parameters:
    - F_y = Member yield strength
    - F_y_min = Minimum allowable probable yield strength"""
    
    if 1.1 * F_y >= F_y_min: R_y = 1.1
    elif 1.1 * F_y < F_y_min: R_y = (F_y_min) / F_y
    
    return R_y

## Beams ##
@handcalc(jupyter_display = False, precision = 2)
def M_b_prob_func(R_y, F_y, Z, R_sh = 1.1):
    """Gives the probable moment resistance of a beam per clause 27.2.2.  The strain hardening factor is 
    equal to 1.1 by default (R_sh = 1.1), but this can be changed as per Annex J.
    
    Parameters:
    - R_y = Probable yield stress factor
    - F_y = Member yield strength
    - Z = Plastic section modulus
    - R_sh = Strain hardening factor (set as 1.1 unless otherwise permitted)"""
    
    M_pb =  Z * F_y
    R_sh
    M_b_prob = R_sh * R_y * M_pb
    
    return M_b_prob

## Columns ##
@handcalc(jupyter_display = False, precision = 2)
def M_prime_rc_capacity(Z, A, F_y, C_f):
    """Gives the column factored flexural resistance projected at the intersection of the beam and column
    centerlines in accordance with clause 27.2.3.3.
    
    Parameters:
    - F_y = Member yield strength
    - Z = Plastic section modulus
    - C_f = Axial force from the gravity loads plus the summation of V_h acting at and above the level under 
    consideration"""
    
    M_pc = Z * F_y
    C_y =  A * F_y
    M_prime_rc = 1.18 * phi * M_pc * (1-C_f/(phi * C_y))
    if M_prime_rc <= M_pc: M_prime_rc = M_prime_rc
    elif M_prime_rc > M_pc: M_prime_rc = M_pc
        
    return M_prime_rc

@handcalc(jupyter_display = False, precision = 2)
def M_prime_rc_demand(M_b_prob, V_h, d_c, x):
    """Gives the probable load acting on the column at the beam-to-column intersection as per clause 27.2.3.3.
    
    Parameters:
    - M_b_prob = 
    - V_h = 
    - d_c = 
    - x = """
    
    M_prime_pc = M_b_prob + V_h*(x+d_c/2)
    
    return M_prime_pc

def M_prime_rc_sum(M_prime_rc_list, M_prime_pc_list):
    """Checks clause 27.2.3.3 by taking the sum of the factored flexural resistance of the column and comparing it to 
    the probable load acting on the column at the beam-to-column intersection.
    
    Parameters:
    - M_prime_rc_list = A list of all factored flexural resistances at the intersection fo the beam and column 
    centerlines
    - M_prime_pc_list = A list of all probable loads acting on the column at the beam-column workpoints"""
    
    Sigma_M_prime_rc = sum(M_prime_rc_list) # sum of Mrc
    # create a string showing the summation of Mrc
    Sigma_M_prime_rc_str = " + ".join(str(x) for x in M_prime_rc_list)
    if len(M_prime_rc_list) > 1: # if the length is greater than 1, show the sum of the list
        Sigma_M_prime_rc_str += "=" + str(Sigma_M_prime_rc)
    display(Math(fr"\sum M'_{{rc}}={Sigma_M_prime_rc_str}")) # display the list
    
    Sigma_M_prime_pc = sum(M_prime_pc_list) # sum of Mpc
    # create a string showing the summation of Mpc
    Sigma_M_prime_pc_str = " + ".join(str(x) for x in M_prime_pc_list)
    if len(M_prime_pc_list) > 1: # if the length is greater than 1, show the sum of the list
        Sigma_M_prime_pc_str += "=" + str(Sigma_M_prime_pc)
    display(Math(fr"\sum M'_{{pc}}={Sigma_M_prime_pc_str}")) # display the list
    
    if Sigma_M_prime_rc >= Sigma_M_prime_pc: result = "ok" # "ok" if the capacity is greater than the demand
    else: result = "not ok" # otherwise not ok
        
    display(Math(fr"\sum M'_{{rc}}={Sigma_M_prime_rc:.1f} \geq \sum M'_{{pc}}={Sigma_M_prime_pc:.1f} " + 
                 fr"\quad \textbf{{{result}}}")) # render math
    
## Joint Panel Zone ##
@handcalc(jupyter_display = False, precision = 2)
def panel_zone_shear_a(d_c, b_c, t_c, d_b, w_prime, F_yc):
    """Gives the horizontal shear resistance of the column joint panel zone according to clause 27.2.4.2 a)
    
    Parameters:
    - d_c = Depth of the column
    - b_c = Width of the column
    - t_c = Thickness of the column flange
    - d_b = Depth of the beam
    - w_prime = Thickness of the column web plus the thickness of the double plates when used
    - F_yc = Column yield strength"""
    
    V_r = 0.55*phi*d_c*w_prime*F_yc*(1+(3*b_c*t_c**2)/(d_c*d_b*w_prime))
    
    if V_r <= 0.66*phi*d_c*w_prime*F_yc: V_r = V_r
    elif V_r > 0.66*phi*d_c*w_prime*F_yc: V_r = 0.66*phi*d_c*w_prime*F_yc
        
    return V_r

@handcalc(jupyter_display = False, precision = 2)
def panel_zone_shear_b(d_c, w_prime, F_yc):
    """Gives the horizontal shear resistance of the column joint panel zone according to clause 27.2.4.2 b)
    
    Parameters:
    - d_c = Depth of the column
    - w_prime = Thickness of the column web plus the thickness of the double plates when used
    - F_yc = Column yield strength"""
    
    V_r = 0.55*phi*d_c*w_prime*F_yc
        
    return V_r

#=== Table 2 Section Class Checks ===#
@handcalc(jupyter_display = False, precision = 1)
def class_check_1_edge(b_el, t, F_y):
    """Gives the class of an element supported along one edge and under flexural compression
    per Table 2, such as: Flanges of I-sections or T-sections about their major axis; Plates 
    projecting from element in compression elements; and Outstanding legs of pairs of angles  in continuous 
    contact with an axis of symmetry in the plane of loading.
    
    Parameters:
    - b_el = Effective element width
    - t = Thickness of the element
    - F_y = Member yield strength
    """
    
    if b_el / t <= 145 / sqrt(F_y): Class = 1
    elif b_el / t <= 170 / sqrt(F_y): Class = 2
    elif b_el / t <= 200 / sqrt(F_y): Class = 3
    elif b_el / t > 200 / sqrt(F_y): Class = 4
    
    return Class

@handcalc(jupyter_display = False, precision = 1)
def class_check_I_web(h, w, F_y, C_f, A):
    """Gives the class of an element supported along two edges and under flexural and axial compression 
    per Table 2, such as Webs of I-sections or T-sections about their major axis.
    
    Parameters:
    - h = Height of the element
    - w = Width of the element
    - F_y = Member yield strength
    - C_f = Factored axial compressive force acting on the the member
    - A = Gross area of the member
    """
    
    C_y = A * F_y
    
    if h / w <= 1100 / sqrt(F_y) * (1 - 0.39 * (C_f / (phi * C_y))): Class = 1
    elif h / w <= 1700 / sqrt(F_y) * (1 - 0.61 * (C_f / (phi * C_y))): Class = 2
    elif h / w <= 1900 / sqrt(F_y) * (1 - 0.65 * (C_f / (phi * C_y))): Class = 3
    elif h / w > 1900 / sqrt(F_y) * (1 - 0.65 * (C_f / (phi * C_y))): Class = 4
    
    return Class

#=== Table 5 ===#
def table_5(d_bt, condition, mm = mm):
    """This function returns the minimum edge distance from Table 5."""
    
    if d_bt <= 1 / 2 * inch:
        sheared, rolled = 26 * mm, 20 * mm
    elif d_bt <= 5 / 8 * inch:
        sheared, rolled = 28 * mm, 22 * mm
    elif d_bt <= 16 * mm:
        sheared, rolled = 28 * mm, 22 * mm
    elif d_bt <= 3 / 4 * inch:
        sheared, rolled = 32 * mm, 25 * mm        
    elif d_bt <= 20 * mm:
        sheared, rolled = 34 * mm, 26 * mm        
    elif d_bt <= 7 / 8 * inch:
        sheared, rolled = 38 * mm, 28 * mm        
    elif d_bt <= 22 * mm:
        sheared, rolled = 38 * mm, 28 * mm        
    elif d_bt <= 24 * mm:
        sheared, rolled = 42 * mm, 30 * mm        
    elif d_bt <= 1 * inch:
        sheared, rolled = 44 * mm, 32 * mm        
    elif d_bt <= 27 * mm:
        sheared, rolled = 48 * mm, 34 * mm
    elif d_bt <= 1 * inch + 1 / 8 * inch:
        sheared, rolled = 51 * mm, 38 * mm
    elif d_bt <= 30 * mm:
        sheared, rolled = 52 * mm, 38 * mm
    elif d_bt <= 1 + 1/4 * inch:
        sheared, rolled = 57 * mm, 41 * mm
    elif d_bt <= 36 * mm:
        sheared, rolled = 64 * mm, 46 * mm
    else:
        sheared, rolled = 1.75 * d_bt, 1.25 * d_bt
    
    if condition == "sheared":
        return sheared
    elif condition == "rolled":
        return rolled