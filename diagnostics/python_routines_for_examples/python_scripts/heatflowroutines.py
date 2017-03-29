#!/usr/bin/python

import sys
import h5py
import numpy as np
import re
from oshunroutines import pullData

def getcompoundxt(directory,quantity,timerange):
    simtime, axis1, axis2, data = pullcompoundData(directory,quantity,0)
    
    dataoveralltime = np.zeros((timerange.shape[0],axis1.shape[0]))
    timeaxis = []
    datamax = []
    it = 0
    for time in timerange:
        simtime, axis1, axis2, data = pullcompoundData(directory,quantity,time)
        dataoveralltime[it,:] = data
        timeaxis.append(simtime)
        datamax.append(max(abs(data)))
        it = it + 1
    
    return timeaxis, axis1, datamax, dataoveralltime
    

def pullcompoundData(directory,quantity,time):
    
    time = int(time)
    if (time < 1e1 and time >= 0):
        timestr = '0000' + str(time);
    elif (time < 1e2 and time >= 1e1):
        timestr = '000' + str(time);
    elif (time < 1e3 and time >= 1e2):
        timestr = '00' + str(time);
    elif (time < 1e4 and time >= 1e3):
        timestr = '0' + str(time);
    else:
        timestr = str(time);
        
    if quantity == 'kappa':
        simtime, axis1, axis2, data_n = pullData(directory,'n',time)
        simtime, axis1, axis2, data_T = pullData(directory,'T',time)
        simtime, axis1, axis2, data_Qx = pullData(directory,'Qx',time)
        simtime, axis1, axis2, data_Z = pullData(directory,'Z',time)
        data_Z = np.ones(len(data_T))
        inputdeck = open(directory + '/inputdeck', 'r') 
        for line in inputdeck: 
            if 'hydrocharge' in line:
                temp = re.findall("\d+\.\d+", line)
                Z0 = float(temp[0])
            if 'density_np' in line:
                temp = re.findall("[\s=]+([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))$", line)
                n0 = float(temp[0])

        data = calcKappa(axis1, data_n, data_T, data_Qx, data_Z*Z0, n0)

    elif quantity == 'kappaEH':
        
        simtime, axis1, axis2, data_Z = pullData(directory,'Z',time)
        simtime, axis1, axis2, data_wt = pullcompoundData(directory,'wt',time)
        
        inputdeck = open(directory + '/inputdeck', 'r') 
        for line in inputdeck: 
            if 'hydrocharge' in line:
                temp = re.findall("\d+\.\d+", line)
                Z0 = float(temp[0])
        
        data_Z = data_Z * Z0
        kappa1, data, kappa3 = getKappaEH(data_wt, data_Z)

    elif quantity == 'nuei':
        
        simtime, axis1, axis2, data_Z = pullData(directory,'Z',time)
        simtime, axis1, axis2, data_n = pullData(directory,'n',time)
        simtime, axis1, axis2, data_T = pullData(directory,'T',time)

        inputdeck = open(directory + '/inputdeck', 'r') 
        for line in inputdeck: 
            if 'hydrocharge' in line:
                temp = re.findall("\d+\.\d+", line)
                Z0 = float(temp[0])
            if 'density_np' in line:
                temp = re.findall("[\s=]+([+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))$", line)
                n0 = float(temp[0])
        
        data_Z = data_Z * Z0
        data_n = data_n * n0

        data = calcnuei(data_n, data_T, data_Z)

    elif quantity == 'ND':

        simtime, axis1, axis2, data_n = pullData(directory,'n',time)
        simtime, axis1, axis2, data_T = pullData(directory,'T',time)

        data_T = np.array(data_T)*511000
        data_n = np.array(data_n)
    
        data = 1.72e9*np.sqrt(data_T**3/(data_n))


    elif quantity == 'wt':

        simtime, axis1, axis2, Bx = pullData(directory,'Bx',time)
        simtime, axis1, axis2, By = pullData(directory,'By',time)
        simtime, axis1, axis2, Bz = pullData(directory,'Bz',time)

        simtime, axis1, axis2, nuei = pullcompoundData(directory,'nuei',time)

        data = np.sqrt(np.square(Bx)+np.square(By)+np.square(Bz))/nuei

    elif quantity == 'wc':

        simtime, axis1, axis2, Bx = pullData(directory,'Bx',time)
        simtime, axis1, axis2, By = pullData(directory,'By',time)
        simtime, axis1, axis2, Bz = pullData(directory,'Bz',time)

        data = np.sqrt(np.square(Bx)+np.square(By)+np.square(Bz))

    return  simtime, axis1, axis2, data

def calcnuei(data_n, data_T, data_Z):

    nuei = np.zeros(len(data_T))
    data_T = np.array(data_T)*511000
    data_n = np.array(data_n)
    

    for ix in range(0,len(data_T)):
        if (data_T[ix] > 10.0*data_Z[ix]**2):
            lnei = max(2.0, 24.0 - 0.5*np.log(data_n[ix]) + np.log(data_T[ix]))
        else:
            lnei = max(2.0, 23.0 - 0.5*np.log(data_n[ix]) + 1.5*np.log(data_T[ix])-np.log(data_Z[ix]))
        
        ND = 1.72e9*np.sqrt(data_T[ix]**3/(data_n[ix]))
        nuei[ix] = np.sqrt(2.0/np.pi)/9.0 * data_Z[ix]*lnei / ND

    return nuei

def calcND(data_n, data_T, data_Z):
    ND = np.zeros(len(data_T))

    data_T = np.array(data_T)*511000
    data_n = np.array(data_n)
    
    ND[ix] = 1.72e9*np.sqrt(data_T**3/(data_n))

    return ND

def calcKappa(x_axis, data_n, data_T, data_Qx, data_Z, n0):


    nuei = calcnuei(data_n*n0, data_T, data_Z)

    data_T=data_T*3
    data_T=np.square(data_T)

    gradT = np.roll(data_T,1) - np.roll(data_T,-1)
    gradT = gradT*0.5/(x_axis[2]-x_axis[1])
    gradT = gradT*data_n/nuei/18

    kappa = data_Qx/gradT

    return kappa
	# data_T = np.array(data_T)*511000
	# data_n = np.array(data_n)
	# data_Qx = np.array(data_Qx)
	# nuei = np.zeros(len(x_axis))

	# for ix in range(0,len(x_axis)):
	# 	if (data_T[ix] > 10.0*Z0**2):
	# 		lnei = max(2.0, 24.0 - 0.5*np.log(data_n[ix]*n0) + np.log(data_T[ix]))
	# 	else:
	# 		lnei = max(2.0, 23.0 - 0.5*np.log(data_n[ix]*n0) + 1.5*np.log(data_T[ix])-np.log(Z0))
		
	# 	ND = 1.72e9*np.sqrt(data_T[ix]**3/(data_n[ix]*n0))
	# 	nuei[ix] = np.sqrt(2.0/np.pi)/9.0 * Z0 * lnei / ND

	# data_T=data_T/511000*3
    


def getKappaEH(wtI,ZI):
    kappapar = np.zeros(len(ZI))
    kappaperp = np.zeros(len(ZI))
    kappaw = np.zeros(len(ZI))


    Z = np.array([1,2,3,4,5,6,7,8,10,12,14,20,30,60,1e6+1])
    g0 = np.array([3.203,4.931,6.115,6.995,7.680,8.231,8.685,9.067,9.673,10.13,10.5,11.23,11.9,12.67,13.58])
    g0p = np.array([6.18,9.3,10.2,9.14,8.6,8.57,8.84,7.93,7.44,7.32,7.08,6.79,6.74,6.36,6.21])
    g1p = np.array([4.66,3.96,3.72,3.6,3.53,3.49,3.49,3.43,3.39,3.37,3.35,3.32,3.3,3.27,3.25])
    c0p = np.array([1.93,1.89,1.66,1.31,1.12,1.04,1.02,0.875,.77,.722,.674,.605,.566,.502,.457])
    c1p = np.array([2.31,3.78,4.76,4.63,4.62,4.83,5.19,4.74,4.63,4.7,4.64,4.65,4.81,4.71,4.81])
    c2p = np.array([5.35,7.78,8.88,8.8,8.8,8.96,9.24,8.84,8.71,8.73,8.65,8.6,8.66,8.52,8.53],)
    g0w = np.array([6.071,15.75,25.65,34.95,43.45,51.12,58.05,64.29,75.04,83.93,91.38,107.8,124.3,145.2,172.7])
    g0pp = np.array([4.01,2.46,1.13,.628,.418,.319,.268,.238,.225,.212,.202,.2,.194,.189,.186])
    g1pp = np.array([2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5])
    c0pp = np.array([.661,.156,.0442,.018,.00963,.00625,.00461,.00371,.003,.00252,.00221,.00185,.00156,.0013,.00108])
    c1pp = np.array([.931,.398,.175,.101,.0702,.0551,.0465,.041,.0354,.0317,.0291,.0256,.0228,.0202,.018])
    c2pp = np.array([2.5,1.71,1.05,.775,.646,.578,.539,.515,.497,.482,.471,.461,.45,.44,.43])

    for x in range(0,len(ZI)):
        
        Zin = ZI[x]
        wt = wtI[x] #/(0.75*np.sqrt(np.pi))
        i=0
        check=1
        
        while check:
            if ((Z[i] < Zin and Z[i+1] > Zin) or (Z[i] == Zin)):
                check = 0
            else:
                i=i+1
           
        intZ = (Zin-Z[i])/(Z[i+1]-Z[i])
        g0a = g0[i] + intZ*(g0[i+1]-g0[i])
        g0pa = g0p[i] + intZ*(g0p[i+1]-g0p[i])
        g1pa = g1p[i] + intZ*(g1p[i+1]-g1p[i])
        c0pa = c0p[i] + intZ*(c0p[i+1]-c0p[i])
        c1pa = c1p[i] + intZ*(c1p[i+1]-c1p[i])
        c2pa = c2p[i] + intZ*(c2p[i+1]-c2p[i])
        g0wa = g0w[i] + intZ*(g0w[i+1]-g0w[i])
        g0ppa = g0pp[i] + intZ*(g0pp[i+1]-g0pp[i])
        g1ppa = g1pp[i] + intZ*(g1pp[i+1]-g1pp[i])
        c0ppa = c0pp[i] + intZ*(c0pp[i+1]-c0pp[i])
        c1ppa = c1pp[i] + intZ*(c1pp[i+1]-c1pp[i])
        c2ppa = c2pp[i] + intZ*(c2pp[i+1]-c2pp[i])

        kappapar[x] = g0a
        kappaperp[x] = (g1pa*wt+g0pa)/(wt**3+c2pa*wt**2+c1pa*wt+c0pa)
        kappaw[x] = wt*(g1ppa*wt+g0ppa)/(wt**3+c2ppa*wt**2+c1ppa*wt+c0ppa)

    return kappapar,kappaperp,kappaw