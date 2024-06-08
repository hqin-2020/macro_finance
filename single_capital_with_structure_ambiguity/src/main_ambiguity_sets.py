import os
import sys
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=200)
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import norm
pd.options.display.float_format = '{:.3g}'.format
sns.set(font_scale = 1.0, rc={"grid.linewidth": 1,'grid.color': '#b0b0b0', 'axes.edgecolor': 'black',"lines.linewidth": 3.0}, style = 'whitegrid')
from datetime import datetime
from tqdm import tqdm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.optimize import fsolve
from scipy.optimize import bisect
now = datetime.now()
plt.rcParams['axes.formatter.useoffset'] = True
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)
sns.set(style="whitegrid", font_scale=1.13, rc={"lines.linewidth": 3.5})


#Parameters
sigma_k1 = 0.0087*np.sqrt(1.4)
sigma_k2 = 0.0038*np.sqrt(1.4)
sigma_z1 = 0
sigma_z2 = 0.055*np.sqrt(1.4)
hat_beta_z = 0.056
phi = 8
delta = 0.01
alpha = 0.092
hat_alpha_z = 0
hat_alpha_k = 0.04
hat_beta_k = 0.04


def extract_eta1s(eta1_contour):
    eta1s = []
    for collection in eta1_contour.collections:
        for path in collection.get_paths():
            points = path.vertices
            eta1s.extend(points)
    return eta1s

def get_alphaz(eta_01, eta_02, eta_11,eta_12):
    alpha_z = (sigma_z1*eta_01+sigma_z2*eta_02 + hat_alpha_z)
    return alpha_z

def get_alphak(eta_01, eta_02, eta_11,eta_12):
    alpha_k = -(sigma_k1*eta_01 + sigma_k2*eta_02 - hat_alpha_k)
    return alpha_k

def get_betaz(eta_01, eta_02, eta_11,eta_12):
    beta_z = hat_beta_z - sigma_z1*eta_11 - sigma_z2*eta_12
    return beta_z

def get_betak(eta_01, eta_02, eta_11,eta_12):
    beta_k = hat_beta_k + sigma_k1*eta_11 + sigma_k2*eta_12
    return beta_k

def find_eta_1_root(q, twoparameter):
    if twoparameter==0:
        rho2 = q**2/((sigma_z1**2 + sigma_z2**2))
        print(f"rho_2 for slope parameters only = {rho2}")
    elif twoparameter==1:
        rho2 = q**2/((sigma_z1**2 + sigma_z2**2)*2) 
        print(f"rho_2 for all four parameters = {rho2}")
        
    def quadratic(eta_11, eta_12, rho2):
         return -hat_beta_z * rho2 + rho2 *(sigma_z1*eta_11 + sigma_z2*eta_12) + (eta_11**2+eta_12**2)/2
    
    def get_contour(quadratic, rho2):
        quadratic_vector = np.vectorize(quadratic)
        X = np.linspace(-5.0, 5.0, 5000)
        Y = np.linspace(-5.0, 5.0, 5000)
        X, Y = np.meshgrid(X, Y)
        W = np.ones((5000,))*rho2
        Z = quadratic_vector(X,Y,W)
        eta1_contour = plt.contour(X, Y, Z, levels=[0],colors='black')
        return eta1_contour
    
    eta1_contour = get_contour(quadratic, rho2)
    eta1s = extract_eta1s(eta1_contour)
    
    eta_01s = []
    eta_02s = []
    eta_11s = []
    eta_12s = []
    
    def eta02_root(eta_02, rho2):
        eta_01 = -(rho2*sigma_z2*eta_02+eta_12*eta_02)/(rho2*sigma_z1+eta_11)
        return (eta_01**2+eta_02**2)/2 + (sigma_z1**2+sigma_z2**2)/2*rho2 - q**2/2 
    

    def eta01_root(eta_02, rho2):
        return -(rho2*sigma_z2*eta_02+eta_12*eta_02)/(rho2*sigma_z1+eta_11)
    
    
    for eta1 in eta1s:
        eta_11 = eta1[0]
        eta_12 = eta1[1]
        eta_02 = bisect(eta02_root, 0, 2,args=(rho2,))
        eta_01 = eta01_root(eta_02, rho2)
        eta_01s.append(eta_01)
        eta_02s.append(eta_02)
        eta_11s.append(eta_11)
        eta_12s.append(eta_12)

    for eta1 in eta1s:
        eta_11 = eta1[0]
        eta_12 = eta1[1]
        eta_02 = bisect(eta02_root, -2, 0,args=(rho2,))
        eta_01 = eta01_root(eta_02, rho2)
        eta_01s.append(eta_01)
        eta_02s.append(eta_02)
        eta_11s.append(eta_11)
        eta_12s.append(eta_12)
        
    min_eta11 = min(eta_11s)
    min_eta12 = min(eta_12s)
    max_eta11 = max(eta_11s)
    max_eta12 = max(eta_12s)
    
    alpha_zs = []
    alpha_ks = []
    beta_zs = []
    beta_ks = []


    for eta_01, eta_02, eta_11,eta_12 in zip(eta_01s, eta_02s, eta_11s, eta_12s):
        alpha_zs.append(get_alphaz(eta_01, eta_02, eta_11,eta_12))
        alpha_ks.append(get_alphak(eta_01, eta_02, eta_11,eta_12))
        beta_zs.append(get_betaz(eta_01, eta_02, eta_11,eta_12))
        beta_ks.append(get_betak(eta_01, eta_02, eta_11,eta_12))
        
    
    return eta1_contour, alpha_zs, alpha_ks, beta_zs, beta_ks

    

solution_q_002_2 = find_eta_1_root(0.2,0)
solution_q_002_4 = find_eta_1_root(0.2,1)


#Plot
plotdir = "./plots/"
fig, ax = plt.subplots(1, 3, figsize=(8, 4))

# First subplot
contour, alpha_zs, alpha_ks, beta_zs, beta_ks = solution_q_002_2
plt.subplot(1, 2, 1)
plt.fill(beta_zs, beta_ks, color='#1874CD', linewidth=2)
plt.ylabel('$\\beta_k$', fontsize='20')
plt.xlabel('$\\beta_1$', fontsize='20')
plt.scatter(beta_zs, beta_ks, color='black', s=4)
plt.scatter(hat_beta_z, hat_beta_k, color='black', s=150)
plt.grid(False)

# Second subplot
plt.subplot(1, 2, 1)
contour, alpha_zs, alpha_ks, beta_zs, beta_ks = solution_q_002_4
plt.fill(beta_zs, beta_ks, color='#983333', linewidth=2)
plt.ylabel('$\\beta_k$', fontsize='20')
plt.xlabel('$\\beta_1$', fontsize='20')
plt.scatter(beta_zs, beta_ks, color='black', s=4)
plt.scatter(hat_beta_z, hat_beta_k, color='black', s=150)
plt.grid(False)

plt.subplot(1, 2, 2)
plt.fill(alpha_zs, alpha_ks, color='#983333', linewidth=2)
plt.ylabel('$\\eta_k$', fontsize='20')
plt.xlabel('$\\phi_1$', fontsize='20')
plt.scatter(alpha_zs, alpha_ks, color='black', s=4)
plt.scatter(hat_alpha_z, hat_alpha_k, color='black', s=150)
plt.grid(False)

xmin = -0.02
xmax = 0.2
ymin = 0.02
ymax = 0.055
ax[0].set_xlim(xmin,xmax)
ax[0].set_ylim(ymin,ymax)
ax[1].set_xlim(xmin,xmax)
ax[1].set_ylim(ymin,ymax)
ax[2].set_xlim(xmin,xmax)
ax[2].set_ylim(ymin,ymax)

plt.tight_layout() 
plt.savefig(plotdir+"figure_4.png")
plt.savefig(plotdir+"figure_4.pdf")