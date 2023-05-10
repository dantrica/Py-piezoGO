#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 16:23:07 2019

@author: David A. Miranda, Ph.D
Edited by: Daniel Andrés Triana Camacho
"""

import numpy as np
from scipy import optimize
import pandas as pd
import matplotlib.pyplot as plt

# == METHOD 2b ==
# Advanced usage, with jacobian
#method_2b  = "leastsq with jacobian"

class semiCirclesEISModel:
    def __init__(self, f, Z):
        self.f = f
        self.k = None
        self.ABC = None
        self.R = None
        self.DeltaR = None
        Z_ = np.array(Z)
        self.real = np.real(Z_)
        self.imag = -np.imag(Z_)
        self.ier = None

    def calc_R(self,xc, yc):
        """ calculate the distance of each 2D points from the center c=(xc, yc) """
        return np.sqrt((self.real-xc)**2 + (self.imag-yc)**2)

    #@countcalls
    def f_2b(self, c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(-xc, -yc, C) """
        A, B, C = c
        x = self.real[self.k]
        y = self.imag[self.k]
        r = x**2 + y**2 + 2*x*A + 2*y*B + C
        return r

    def getFittedCircumference(self, x = None):
        A, B, C = self.ABC
        R = self.R
        #print 'x', x
        if type(x) == type(None):
            x = self.real
        k = (R**2 - (x + A)**2) >= 0
        y = np.array(x.size * [0.0])
        y1 = np.array(x.size * [0.0])
        y2 = np.array(x.size * [0.0])
        y1[k] = -B + np.sqrt(R**2 - (x[k] + A)**2)
        y2[k] = -B - np.sqrt(R**2 - (x[k] + A)**2)
        r1 = np.sum((self.imag - y1)**2)
        r2 = np.sum((self.imag - y2)**2)
        y = y1
        #print r1, r2
        if r2 < r1:
            y = y2
        return x, y

    #@countcalls
    def Df_2b(self, c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        A, B, C     = c
        xc, yc = -A, -B
        x = self.real[self.k]
        y = self.imag[self.k]
        df2b_dc    = np.empty((len(c), x.size))

        Ri = np.sqrt(A**2 + B**2 - C) #self.calc_R(xc, yc)
        df2b_dc[ 0] = (xc - x)/Ri                   # dR/dxc
        df2b_dc[ 1] = (yc - y)/Ri                   # dR/dyc
        df2b_dc[ 2] = (xc**2 + yc**2 - Ri)/Ri               # dR/dyc
        df2b_dc       = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

        return df2b_dc

    def optCircleParameters(self, real_min=0.0, real_max=np.Inf, maxNumDataToFit=np.Inf):
        # coordinates of the barycenter
        k = ( self.real >= real_min ) & ( self.real <= real_max )
        k = np.arange(self.real.size)[k]
        k2 = np.linspace(0, k.size-1, int(np.min([maxNumDataToFit, k.size])))
        k2 = np.int16(k2)
        self.k = k[k2]
        x_m = np.mean(self.real[k])
        R0  = (self.real[k].max() - self.real[k].min())/2.0
        y_m = self.imag[k].max() - R0
        C0  =  x_m**2 + y_m**2 - R0**2

        ABC = -x_m, -y_m, C0
        ABC_opt, ier = optimize.leastsq(self.f_2b, ABC)
        A, B, C = ABC_opt
        Ri_2b        = np.sqrt(A**2 + B**2 - C) #self.calc_R(xc_2b, yc_2b)
        
        self.ier = ier
        self.ABC = ABC_opt
        self.R   = Ri_2b
        self.DeltaR = 2*np.sqrt(np.abs(B**2 - Ri_2b**2))
        return ABC_opt, Ri_2b

    def optLinearParameters(self, real_min=0.0, real_max=np.Inf):
        s = ( self.real >= real_min ) & ( self.real <= real_max )
        x = self.real[s]
        y = self.imag[s]
        f = self.f[s]
        p = np.polyfit(x, y, deg=1)
        return
    def search_min_max(self):
        l = np.int16(np.linspace(0, self.real.size-1, 30))
        der = np.gradient(self.imag[l])/np.gradient(self.real[l])
        real_min = []
        real_max = []
        for j1 in np.arange(len(der)-2):
            if (der[j1]*der[j1+1]<0) & (der[j1]*der[j1+2]<0):
                if der[j1]<0:
                    s = der == der[j1]
                    real_min.append(self.real[l][s][0])
                elif der[j1]>0:
                    s = der == der[j1]
                    real_max.append(self.real[l][s][0])
        if len(real_min) == 0:
            real_min.append(self.real[l].max())
        return real_min, real_max, l

    def optTwoDispersions(self, niter = 1000):
        x = self.real
        y = self.imag
        ABC1, R1 = self.optCircleParameters(real_max = x.mean())
        ABC2, R2 = self.optCircleParameters(real_min = x.mean())
        p = np.append(ABC1, ABC2)

        def y_fit(x, ABC):
            A, B, C = ABC
            R = np.sqrt(A**2 + B**2 - C)
            return -B + np.sqrt(R**2 - (x + A)**2)

        def y2disp(x, p):
            p1 = p[:3]
            p2 = p[3:]
            return y_fit(x, p1) + y_fit(x, p2)

        def r_y(p):
            A1, B1, C1, A2, B2, C2 = p
            if A1**2 + B1**2 - C1 < 0 or A2**2 + B2**2 - C2 < 0:
                return 100.0
#            y1 = y_fit(x, [A1, B1, C1])
#            y2 = y_fit(x, [A2, B2, C2])
#            if (y1 < 0).any() or (y2 < 0).any():
#                return 110.0
            return np.sqrt(np.sum( (y - y2disp(x, p))**2 ))/y.size

        p_opt = optimize.basinhopping(r_y, p, niter=niter).x

        self.plotNyquist(fitted=False)
        plt.plot(x, y2disp(x, p_opt), label='modelo de suma de dos semicirculos')
        plt.plot(x, y_fit(x, p_opt[3:]), label='una de las dispersiones')
        plt.plot(x, y_fit(x, p_opt[:3]), label='otra de las dispersiones')
        plt.legend()

    def plotEIS(self, real_marker = 'k', imag_marker = 'k'):
        plt.subplot(211)
        plt.semilogx(self.f, self.real, real_marker)
        plt.ylabel(r'$real\{Z\}$ / $\Omega$')
        plt.xlabel(r'Frequency / $Hz$')

        plt.subplot(212)
        plt.semilogx(self.f, self.imag, real_marker)
        plt.ylabel(r'$-imag\{Z\}$ / $\Omega$')
        plt.xlabel(r'Frequency / $Hz$')

    def plotNyquist(self, data_marker = 'ko', fit_marker = 'k:',
                    data_label = 'Experimental Data',
                    data=True,fitted=True,
                    fit_label='Fitted circumference', show_title = False, data_alpha=0.7, show_legends = True):
        swOpt = False
        if type(self.ABC) != type(np.array([])):
            swOpt = True
        if type(self.k) != type(np.array([])):
            swOpt = True
        if swOpt:
            self.optCircleParameters()

        x = self.real
        y = self.imag
        k = self.k
        A, B, C = self.ABC
        R = np.sqrt(A**2 + B**2 - C)

        q = np.linspace(0,2*np.pi,1000)
        x_ = R*np.cos(q) - A
        y_ = R*np.sin(q) - B

        k = np.int16(np.linspace(0, x.size-1, 30))

        if data:
            plt.plot(x[k], y[k], data_marker, label=data_label, alpha=data_alpha)
        if fitted:
            plt.plot(x_, y_, fit_marker, lw=1, label=fit_label)

        if show_legends:
            plt.legend()

        plt.xlabel(r'$real\{Z\}$ / $\Omega$')
        plt.ylabel(r'$-imag\{Z\}$ / $\Omega$')

        if fitted and show_title:
            plt.title('Center: (%0.2f,%0.2f). R = %0.2f' % (-A, -B, R))

def parseCompexNumber(data):
    if type(data) == type(pd.Series([])):
        return data.apply(strToComplex)
    if type(data) == type('str'):
        return strToComplex(data)
    return data

def strToComplex(cStr):
    try:
        realAndImag = []
        by = ''
        if cStr.find('+') > 0:
            by = '+'
        if cStr.find('-') > 0:
            by = '-'
        if by != '':
            realAndImag = cStr.split(by)
            realAndImag[1] = (by if by == '-' else '') + realAndImag[1][:-1]
        elif cStr.find('i'):
            realAndImag = ['0',by+cStr[:-1]]
        else:
            realAndImag = [by+cStr[:-1],'0']
        return np.float(realAndImag[0]) + 1j*np.float(realAndImag[1])
    except:
        return cStr
