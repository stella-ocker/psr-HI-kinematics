import astropy
import numpy as np
import pandas as pd
import csv
from kinematic_functions import *

name = input('Source name: ')
ls = input('Galactic longitude (deg): ')
bs = input('Galactic latitude (deg): ')
vls = input('Lower velocity bound (km/s), input nan if none, T if tangent: ')
vus = input('Upper velocity bound (km/s), input nan if none, T if tangent: ')
verrs = input('Velocity error (km/s): ')
plots = input('Make plots? (y/n): ') # currently only works when both distance bounds available

have_lower_bound = True # initialize
have_upper_bound = True # initialize

l = np.deg2rad(float(ls))
b = np.deg2rad(float(bs))

if vls != 'nan' and vls!='T':
    vl = float(vls)
if vus !='nan' and vus!='T':
    vu = float(vus)
if vls == 'nan':
    have_lower_bound = False
if vus == 'nan':
    have_upper_bound = False
verr = float(verrs)

distances = np.linspace(0, 20, 100000)

radius_reid = reid_galactocentric_radius(l)
omega_reid = angular_velocities_translation(calculate_universal_rotation_curve(radius_reid), radius_reid)
velocities_reid = get_radial_velocities(omega_reid, l, b)
lower_intersection_reid = (0,0)
upper_intersection_reid = (0,0)

d_l = 0
d_u = 0
dl_err = 0
du_err = 0
        
if have_lower_bound:
    if vls == 'T':
        print('Lower velocity at tangent point')
        d_l = R0 * np.cos(l)
        ind_tp = np.where(distances>=d_l)[0][0]
        low_vel = velocities_reid[ind_tp]-verr
        lower_intersection = calculate_intersection(low_vel,velocities_reid)
        dl_err = abs(d_l - lower_intersection[0])
        print('Lower distance = ',f"{d_l:.4f}",' +/- ',f"{dl_err:.4f}",' kpc')

    elif np.rad2deg(l) >= 100 and np.rad2deg(l) <= 150:
        print('Applying Perseus correction to lower velocity')
        lower_intersection_reid = calculate_intersection(vl/1.6, velocities_reid)
        lower_intersection_err = calculate_intersection((vl/1.6)-verr, velocities_reid)
        d_l = lower_intersection_reid[0]
        dl_err = abs(lower_intersection_reid[0] - lower_intersection_err[0])
        print('Lower distance = ',f"{d_l:.4f}",' +/- ',f"{dl_err:.4f}",' kpc')
            
    else:
        print('Lower velocity bound available')
        lower_intersection_reid = calculate_intersection(vl, velocities_reid)
        lower_intersection_err = calculate_intersection(vl-verr,velocities_reid)
        d_l = lower_intersection_reid[0]
        dl_err = abs(lower_intersection_reid[0] - lower_intersection_err[0])
        print('Lower distance = ',f"{d_l:.4f}",' +/- ',f"{dl_err:.4f}",' kpc')

if have_upper_bound:
    if vus == 'T':
        print('Upper velocity at tangent point')
        d_u = R0 * np.cos(l)
        ind_tp = np.where(distances>=d_u)[0][0]
        low_vel = velocities_reid[ind_tp]-verr
        lower_intersection = calculate_intersection(low_vel,velocities_reid)
        du_err = abs(d_u - lower_intersection[0])
        print('Upper distance = ',f"{d_u:.4f}",' +/- ',f"{du_err:.4f}",' kpc')
            
    if np.rad2deg(l) >= 100 and np.rad2deg(l) <= 150:
        print('Applying Perseus correction to upper velocity')
        upper_intersection_reid = calculate_intersection(vu/1.6, velocities_reid)
        upper_intersection_err = calculate_intersection((vu/1.6)-verr, velocities_reid)
        d_u = upper_intersection_reid[0]
        du_err = abs(upper_intersection_reid[0] - upper_intersection_err[0])
        print('Upper distance = ',f"{d_u:.4f}",' +/- ',f"{du_err:.4f}",' kpc')
            
    else:
        print('Upper velocity bound available')
        upper_intersection_reid = calculate_intersection(vu, velocities_reid)
        upper_intersection_err = calculate_intersection(vu-verr,velocities_reid)
        d_u = upper_intersection_reid[0]
        du_err = abs(upper_intersection_reid[0] - upper_intersection_err[0])
        print('Upper distance = ',f"{d_u:.4f}",' +/- ',f"{du_err:.4f}",' kpc')

if plots and have_upper_bound and have_lower_bound and vus != 'T' and vls != 'T':

    plt.figure(figsize=(5,3))
    plt.plot(velocities_reid,distances,lw=2,color='dodgerblue',label='R19 Rotation Curve')
    plt.vlines(vl,0,20,color='tab:orange',label=r'$V_{\rm HI}$ Bounds')
    plt.vlines(vu,0,20,color='tab:orange')
    plt.ylim(20,0)
    plt.legend(loc='lower right',framealpha=1)
    plt.ylabel('Distance (kpc)',fontsize=11)
    plt.xlabel('Radial Velocity (km s$^{-1}$)',fontsize=11)
    plt.title(name)
    plt.tight_layout()
    plt.tick_params(which='major',direction='in',top=True,right=True)
    plt.show()