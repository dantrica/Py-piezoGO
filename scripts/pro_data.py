# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 23:33:44 2019

@author: Silvi

20210128 update python 2 code to python 3 code 
We change dictionary.has_key('key'): for 'key' in dictionary:
20210129 We need to change the attribute join from python 2.0 instead attribute in python 3.0

    name = '' # Name of experiment
    date = strftime('%c')
    F = [] #Force in N
    ds = [] #deformation in m
    mec_rate = [] #rate change of the force in time
    f_charge = []
    path_m = [] #file name of mechanical proofs
    path_is = [] #file name of impedance proofs
    path_ocp = [] #file name of OCP proofs
    freq = [] #frequency in Hz
    Z_re = [] #real part of impedance in Ohm
    Z_im = [] #imaginary part of impedance in Ohm
    Z = [] #magnitude of impedance in Ohm
    phase = [] #phase of impedance
    V_ocp = [] #open circuit potential
    t_ocp = [] #time ocp measurement
    t_force = [] #time of compression test
    scvPhi = []  # External potential for SCV experiment
    t = [] # Time step
    V = [] # Potential difference
    I = [] # Current
"""

import numpy as np
import matplotlib.pyplot as plt
from time import strftime
from scipy import stats
from scipy.interpolate import interp1d
from lmfit import Model
import string
from scipy.integrate import simps

class pro_data:
    def __init__(self):
        self.name = ''
        self.date = ''
        self.F = []
        self.ds = []
        self.mec_rate = []
        self.f_charge = []
        self.path_m = []
        self.path_is = []
        self.path_ocp = []
        self.set_date()
        self.freq = []
        self.Z_re = []
        self.Z_im = []
        self.Z = []
        self.phase = []
        self.V_ocp = []
        self.t_ocp = []
        self.t_force = []
        self.C_is = []
###################
        self.scvPhi = []
        self.t = []
        self.V = []
        self.I = []
#########################
    def set_date(self):
        self.date = strftime('%c')
    def select_data(self, json):
        if type(json) != type({'1':''}):
            print('Error: the input object must be a dictionary')
        else:
            self.set_date()
            if 'name' in json:
            # if json.has_key('name'):
                self.name = json['name']
            if 'compresion' in json:
            # if json.has_key('compresion'):
                for mec_data in json['compresion']:
                    self.add_mec_data(mec_data)
            if 'is' in json:
                for imp_data in json['is']:
                    self.add_imp_data(imp_data)
            if 'ocp' in json:
            # if json.has_key('ocp'):
                for ocp_data in json['ocp']:
                    self.add_ocp_data(ocp_data)
##############################################3
            if 'scv' in json:
            # if json.has_key('scv'):
                for scvData in json['scv']:
                    self.addSCVdata(scvData)
#################################################3
    def add_mec_data(self, mec_data):
        if 'path' in mec_data and 'rate' in mec_data:
        # if mec_data.has_key('path') and mec_data.has_key('rate'):
            self.path_m.append(mec_data['path'])
            data_c = np.loadtxt(fname=mec_data['path'], skiprows=1)
            self.mec_rate.append(mec_data['rate'])
            self.F.append(data_c[:,0])
            self.ds.append(data_c[:,1])
            self.t_force.append(data_c[:,2])
        else:
            print('Error: path or rate keys were omited')
    def add_ocp_data(self, ocp_data):
        if 'path' in ocp_data and 'rate' in ocp_data:
        #if ocp_data.has_key('path') and ocp_data.has_key('rate'):
            self.path_ocp.append(ocp_data['path'])
#            data_o = np.loadtxt(fname=ocp_data['path'], skiprows=1, usecols=(0,1,2,3,4))
            data_o = np.loadtxt(fname=ocp_data['path'], skiprows=1)
            self.t_ocp.append(data_o[:,0])
            self.V_ocp.append(data_o[:,1])
        else:
            print('Error: path or rate keys were omited')
    def add_imp_data(self, imp_data):
        #if 'path' in imp_data and 'force' in imp_data:
        if 'path' in imp_data:
        #if imp_data.has_key('path') and imp_data.has_key('force'):
            self.path_is.append(imp_data['path'])
            data_is = np.loadtxt(fname=imp_data['path'], skiprows=1)
            try:
                self.f_charge.append(imp_data['force'])
            except:
                pass
            self.freq.append(data_is[:,1])
            self.Z_re.append(data_is[:,2])
            self.Z_im.append(data_is[:,3])
            self.Z.append(data_is[:,4])
            self.phase.append(data_is[:,5])
            self.C_is.append(1/(2*np.pi*self.freq[-1]*(self.Z_re[-1]-1j*self.Z_im[-1])))
        else:
            print('Error: path or force keys were omited')
##########################
    def addSCVdata(self, scvData):
        if 'path' in scvData and 'phi' in scvData:
        #if scvData.has_key('path') and scvData.has_key('phi'):
            data = np.loadtxt(fname=scvData['path'], skiprows=1, usecols=(0,1,2,3,4,5,6))
            self.scvPhi.append(scvData['phi'])
            self.t.append(data[:,1])
            self.I.append(data[:,4])
            self.V.append(data[:,5])
        else:
            print('Error: path or phi keys were omited')
###########################
    def get_force(self):
        rate = np.array([])
        tf = np.array([])
        strain = np.array([])
        force = np.array([])
        ct_force = np.array([])
        ct_ocp = np.array([])
        slope = 0
        intercept = 0
        for j in range(len(self.F)):
            ct_force = np.array(self.t_force)[j] - np.array(self.t_force)[j][0]
            iforce = interp1d(ct_force, self.F[j])
            strain = np.append(strain, self.ds[j])
            rate = np.append(rate, self.mec_rate)
            #condicion para usar solo los datos de fuerza
            if len(self.t_ocp) == 0:
                ct_ocp = ct_force
            else:
                ct_ocp = np.array(self.t_ocp)[j] - np.array(self.t_ocp)[j][0]
            #condicion para realizar la interpolacion con el vector de tiempos mas pequeño
            if ct_force[-1] < ct_ocp[-1]:
                force = np.append(force, iforce(ct_force))
                tf = np.append(tf, ct_force)
                slope, intercept, r_value, p_value, std_err = stats.linregress(ct_force, abs(iforce(ct_force)))
            else:
                force = np.append(force, iforce(ct_ocp))
                tf = np.append(tf, ct_ocp)
                slope, intercept, r_value, p_value, std_err = stats.linregress(ct_ocp, abs(iforce(ct_ocp)))
        return slope, intercept, rate, strain, tf, abs(force)
    def get_ocp(self):
        tocp = np.array([])
        Vocp = np.array([])
        ct_ocp = np.array([])
        ct_force = np.array([])
        for j in range(len(self.V_ocp)):
            ct_ocp = np.array(self.t_ocp)[j] - np.array(self.t_ocp)[j][0]
            iocp = interp1d(ct_ocp, self.V_ocp[j])#[::-1])#para p11
            #condicion para usar solo los datos de ocp
            if len(self.t_force) == 0:
                ct_force = ct_ocp
            else:
                ct_force = np.array(self.t_force)[j] - np.array(self.t_force)[j][0]
            #condicion para realizar la interpolacion con el vector de tiempos mas pequeño
            if ct_ocp[-1] < ct_force[-1]:
                Vocp = np.append(Vocp, iocp(ct_ocp))
                tocp = np.append(tocp, ct_ocp)
            else:
                Vocp = np.append(Vocp, iocp(ct_force))
                tocp = np.append(tocp, ct_force)
        return tocp, Vocp
    def get_Z(self):
        charge = self.f_charge
        f = self.freq
        Zm = self.Z
        Zr = self.Z_re
        Zi = self.Z_im
        return charge, f, Zm, Zr, Zi
    
    def get_C_is(self):
        f = self.freq
        C = self.C_is
        return f, C
############################
    def getC_mu_ti(self, ti):
        phi = np.array([])
        V_t = np.array([])
        Q = np.array([])
        # doc = open('Capacitancia ti = '+str(1e3*ti)+'ms.txt','w+')
        # data = string.join(['index','Potential','Capacitance'],sep='\t')
        # doc.write(data+'\n')
        for k in range(len(self.I)):
            fV = interp1d(self.t[k], self.V[k])
            index_ti = self.t[k] <= ti # ti: is time integration
            q = simps(self.I[k][index_ti], self.t[k][index_ti])
            phi  = np.append(phi, self.scvPhi[k])
            V_t = np.append(V_t, fV(ti))
            Q = np.append(Q, q)
        C_mu_ti = np.zeros(Q.shape, np.float)
        C_mu_ti[0:-1] = np.gradient(Q)/np.gradient(V_t)
        # for j in range(len(C_mu_ti)):
            # data = string.join([str(j+1),str(V_t[j]).replace(".",","),str(1e6*C_mu_ti[j]).replace(".",",")],sep='\t')
            # doc.write(data+'\n')
        return V_t, phi, C_mu_ti, Q
##################################
    def plot_force(self):
        slope, intercept, rate, strain, tf, force = self.get_force()
        plt.figure(1)
        label_r = r'$m_F = %0.6f [KN/m^2], b_F = %0.5f [kN]$' % (1000*slope, intercept)
        plt.plot(-np.flip(strain), np.flip(force), '.')
        #plt.plot(tf, slope*tf+intercept, linestyle = '-', label = label_r)
        plt.xlabel(r'Strain [mm]', fontsize=18)#(r'Axial Displacement [mm]')
        plt.ylabel(r'Force [kN]', fontsize=18)
        #plt.title()
        plt.legend()
        plt.show()
    def plot_ocp(self):
        plt.figure(2)
        for j in range(len(self.path_ocp)):
            plt.plot(self.t_ocp[j], abs(self.V_ocp[j])*1e3, '.', label = self.path_ocp[j][0:18])# + ', ' + str(self.f_charge[j]) + str(' kN'))        
        plt.xlabel(r'Time [s]')
        plt.ylabel(r'Voltage [mV]')
        #plt.title()
        plt.legend()
        plt.show()
#     def plot_ocp_force(self):
#         tocp, Vocp = self.get_ocp()
#         slope, intercept, rate, strain, tf, force = self.get_force()
#         m_VF, b_VF, r_value, p_value, std_err = stats.linregress(force, Vocp)
#         label_r = r'$m = %0.6f [mV/kN], b = %0.5f [mV], r^2=%0.5f$' % (1e3*m_VF, 1e3*b_VF, r_value)
#         label_data = 'Sample 15' #POR MEJORAR
#         plt.figure(3)
# #        plt.subplot(326)
# #        plt.cla()
#         plt.plot(force, Vocp*1e3, '.', label=label_data)
#         plt.plot(force, (m_VF*force+b_VF)*1e3, '-k', label=label_r)
#         plt.xlim(0,2.0)
#         plt.ylim(-0.164*1e3, 0.026*1e3)
#         plt.xlabel(r'Force [kN]', fontsize=18)
#         plt.ylabel(r'OCP [mV]', fontsize=18)
# #        plt.title(r'Reference (curing)    (a)')
# #        plt.title(r'Reference (curing+electric field) (b)')
# #        plt.title(r'442 ppm AuNP (curing) (c)')
# #        plt.title(r'442 ppm AuNP (curing+electric field) (d)')
# #        plt.title(r'658 ppm AuNP (curing)    (e)')
# #        plt.title(r'658 ppm AuNP (curing+electric field) (f)')
#         plt.legend(loc=3, fontsize=12)
#         plt.show()
    def plot_Z(self, color = '.k'):
        charge, f, Zm, Zre, Zim = self.get_Z()
        # plt.figure()
        # plt.xscale('log')
        # plt.xlabel(r'Frequency [Hz]', fontsize=18)#plt.xlabel(r'Frequency [Hz]')
        # plt.ylabel(r'Z [$k\Omega$]', fontsize=18)
        # plt.title(r'Reference (curing)    (a)')
        # plt.title(r'Reference (curing+electric field) (b)')
        # plt.title(r'442 ppm AuNP (curing) (c)')
        # plt.title(r'442 ppm AuNP (curing+electric field) (d)')
        # plt.title(r'658 ppm AuNP (curing)    (e)')
        # plt.title(r'658 ppm AuNP (curing+electric field) (f)')
        # plt.legend()
        plt.figure()
#        plt.subplot(321)
        for s in range(len(self.freq)):
            plt.plot(Zre[s], Zim[s], '.', label = self.path_is[s][0:19])#[0:16])
        #plt.xlim(0, 12e3*1e-3)
        #plt.ylim(0, 12e3*1e-3)
        plt.xlabel(r"real{Z} [$k\Omega$]", fontsize=18)
        plt.ylabel(r'imag{Z} [$k\Omega$]', fontsize=18)
#        plt.title(r'Reference (curing)    (a)')
#        plt.title(r'Reference (curing+electric field) (b)')
#        plt.title(r'442 ppm AuNP (curing) (c)')
#        plt.title(r'442 ppm AuNP (curing+electric field) (d)')
#        plt.title(r'658 ppm AuNP (curing)    (e)')
#        plt.title(r'658 ppm AuNP (curing+electric field) (f)')
        plt.legend()
        # plt.show()
    def plot_C_is(self):
        f, C = self.get_C_is()
        plt.figure()
        for i in range(len(self.freq)):
            label = self.path_is[i]
            plt.subplot(2,1,1)
            plt.plot(f[i], C[i].real*1e6, '.', label=label)
            plt.xscale('log')
            plt.xlabel(r'Frequency [Hz]', fontsize=18)
            plt.ylabel(r'real{C} [$\mu F$]', fontsize=18)
            plt.legend()
            plt.subplot(2,1,2)
            plt.plot(f[i], C[i].imag*1e6, '.', label=label)
            plt.xscale('log')
            plt.xlabel(r'Frequency [Hz]', fontsize=18)
            plt.ylabel(r'imag{C} [$\mu F$]', fontsize=18)
            plt.legend()
            plt.show()
        plt.figure()
        for i in range(len(self.freq)):
            label = self.path_is[i]
            plt.plot(C[i].real*1e6, C[i].imag*1e6, '.', label=label)
            plt.xlabel(r'real{C} [$\mu F$]', fontsize=18)
            plt.ylabel(r'imag{C} [$\mu F$]', fontsize=18)
            plt.legend()
            plt.show()
#####################################################
    def plotC_mu_ti(self, ti, color = '.b', n = 0):
        V_t, phi, C_mu_ti, Q = self.getC_mu_ti(ti)
        #nombre = self.setData()
        #n = 1 name of solution in title and fsample is show in label
        if n == 0:
            name_in_label = str(self.name)
            name_in_title = '$t_{integration}$ = '+str(1e3*ti)+' [ms]'
        else:
            name_in_label = '$t_{integration}$ = '+str(1e3*ti)+' [ms]'
            name_in_title = str(self.name)
        plt.plot(V_t, 1e6*C_mu_ti, color, label = name_in_label)
        plt.xlabel(r'$\phi [V]$')
        plt.ylabel(r'$C [\mu F]$')
        plt.title(name_in_title)
        plt.legend()
        #plt.show()
########################################################