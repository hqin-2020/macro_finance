import mfr.modelSoln as m
import numpy as np
import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description="parameter settings")
parser.add_argument("--nW",type=int,default=30)
parser.add_argument("--nZ",type=int,default=30)
parser.add_argument("--dt",type=float,default=0.1)
parser.add_argument("--a_e",type=float,default=0.14)
parser.add_argument("--a_h",type=float,default=0.135)
parser.add_argument("--rho_e",type=float,default=1.0)
parser.add_argument("--rho_h",type=float,default=1.0)
parser.add_argument("--gamma_e",type=float,default=1.0)
parser.add_argument("--gamma_h",type=float,default=1.0)
parser.add_argument("--chiUnderline",type=float,default=1.0)
parser.add_argument("--action_name",type=str,default='')
parser.add_argument("--shock_expo",type=str,default='')
parser.add_argument("--delta_e",type=float,default=0.01)
parser.add_argument("--delta_h",type=float,default=0.01)
parser.add_argument("--lambda_d",type=float,default=0.02)
parser.add_argument("--nu",type=float,default=0.01)
args = parser.parse_args()

params = m.paramsDefault.copy()

muz2 = 0.0000063030303030303026
## Dimensionality params
params['nDims']             = 2
params['nShocks']           = 2

## Grid parameters 
params['numSds']            = 5
params['uselogW']           = 1
params['nWealth']           = args.nW
params['nZ']                = args.nZ
params['nV']                = 0
params['nVtilde']           = 0

# Domain params
params['Z_bar']             = 0.0
params['V_bar']             = 1.0
params['Vtilde_bar']        = 0.0
params['sigma_K_norm']      = 1*np.sqrt(12)*np.sqrt(muz2)
params['sigma_Z_norm']      = 2.30608/0.4027386142660167*np.sqrt(12)*np.sqrt(muz2)
params['sigma_V_norm']      = 0.0
params['sigma_Vtilde_norm'] = 0.0
params['wMin']              = 0.001
params['wMax']              = 0.999

## Economic params
params['nu_newborn']        = args.nu
params['lambda_d']          = args.lambda_d
params['lambda_Z']          = 0.056
params['lambda_V']          = 0.0
params['lambda_Vtilde']     = 0.0
params['delta_e']           = args.delta_e
params['delta_h']           = args.delta_h
params['a_e']               = args.a_e
params['a_h']               = args.a_h
params['rho_e']             = args.rho_e
params['rho_h']             = args.rho_h
params['phi']               = 8.0
params['gamma_e']           = args.gamma_e
params['gamma_h']           = args.gamma_h
params['equityIss']         = 2
params['chiUnderline']      = args.chiUnderline
params['alpha_K']           = 0.04

## Alogirthm behavior and results savings params
params['method']            = 2
params['dt']                = args.dt
params['dtInner']           = args.dt

params['tol']               = 1e-5
params['innerTol']          = 1e-5

params['verbatim']          = -1
params['maxIters']          = 500000
params['maxItersInner']     = 500000
params['iparm_2']           = 28
params['iparm_3']           = 0
params['iparm_28']          = 0
params['iparm_31']          = 0
params['overwrite']         = 'Yes'
params['exportFreq']        = 1000000
params['CGscale']           = 1.0
params['hhCap']             = 1
params['precondFreq']       = -1

## Shock correlation params
if args.shock_expo == 'lower_triangular':
    params['cov11']             = 1.0
    params['cov12']             = 0.0
    params['cov13']             = 0.0
    params['cov14']             = 0.0
    params['cov21']             = 0.4027386142660167
    params['cov22']             = 0.9153150324227657
    params['cov23']             = 0.0
    params['cov24']             = 0.0
    params['cov31']             = 0.0
    params['cov32']             = 0.0
    params['cov33']             = 0.0
    params['cov34']             = 0.0
    params['cov41']             = 0.0
    params['cov42']             = 0.0
    params['cov43']             = 0.0
    params['cov44']             = 0.0
elif args.shock_expo == 'upper_triangular':
    params['cov11']             = 0.9153150324227657
    params['cov12']             = 0.4027386142660167
    params['cov13']             = 0.0
    params['cov14']             = 0.0
    params['cov21']             = 0.0
    params['cov22']             = 1.0
    params['cov23']             = 0.0
    params['cov24']             = 0.0
    params['cov31']             = 0.0
    params['cov32']             = 0.0
    params['cov33']             = 0.0
    params['cov34']             = 0.0
    params['cov41']             = 0.0
    params['cov42']             = 0.0
    params['cov43']             = 0.0
    params['cov44']             = 0.0

rho_e = "{:0.3f}".format(params['rho_e'])
rho_h = "{:0.3f}".format(params['rho_h'])
gamma_e = "{:0.3f}".format(params['gamma_e'])
gamma_h = "{:0.3f}".format(params['gamma_h'])
a_e = "{:0.3f}".format(params['a_e'])
a_h = "{:0.3f}".format(params['a_h'])
chiUnderline = "{:0.3f}".format(params['chiUnderline'])
dt = "{:0.3f}".format(params['dt'])
delta_e = "{:0.3f}".format(params['delta_e'])
delta_h = "{:0.3f}".format(params['delta_h'])
lambda_d = "{:0.3f}".format(params['lambda_d'])
nu = "{:0.3f}".format(params['nu_newborn'])

folder_name = ('output/' + args.action_name + '/'+args.shock_expo+'/dt_'+str(args.dt)+'/nW_'+str(args.nW)+'_nZ_'+str(args.nZ)+'/chiUnderline_' + chiUnderline + '/a_e_' + a_e + '_a_h_' + a_h  + '/gamma_e_' + gamma_e + '_gamma_h_' + gamma_h + '/rho_e_' + rho_e + '_rho_h_' + rho_h + '/delta_e_' + delta_e + '_delta_h_' + delta_h + '/lambda_d_' + lambda_d + '_nu_' + nu)

params['folderName']        = folder_name
params['preLoad']           = folder_name
Model = m.Model(params)

Model.solve()
Model.printInfo() 
Model.printParams() 
Model.computeStatDent()
Model.dumpData()


