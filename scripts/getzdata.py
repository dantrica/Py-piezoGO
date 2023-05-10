# -*- coding: utf-8 -*-
import os
import numpy as np
class Datos(object):
    def __init__(self,dirname=[], filename=[],
    Units=['[Z in Ohm], [frequencies in Hz], [Phase in degrees]'],
    Index=[],freqs=[],realZ=[],imagZ=[],magnZ=[],phasZ=[],Times=[],Data=[]):
        self.dirname = dirname 
        self.filename = filename
        self.Units = Units
        self.Index = Index
        self.freqs = freqs 
        self.realZ = realZ
        self.imagZ = imagZ
        self.magnZ = magnZ
        self.phasZ = phasZ
        self.Times = Times
        self.Data  = Data
        self.p_Colefit = []
        self.p_multiColefit = []
        self.p_GEMTIPsphfit = []
    #--- Función para llenar datos en una clase ---
    def rellenar(self,DirN,FileN, Units, Z):
        """Llena los datos en a partir de
        DirN:  Directorio base (str)
        FileN: Nombre del archivo (str)
        Units: Cabecera de datos (str)
        Z:     Valores de impedancia (numpy.array)
        """            
        self.dirname.insert(0,DirN)
        self.filename.insert(0,FileN)
        self.Units.insert(0,Units)
        self.Index.insert(0,Z[0])
        self.freqs.insert(0,(Z[1].astype(np.float)))
        self.realZ.insert(0,(Z[2].astype(np.float)))
        self.imagZ.insert(0,(Z[3].astype(np.float)))
        self.magnZ.insert(0,(Z[4].astype(np.float)))
        self.phasZ.insert(0,Z[5])
        self.Times.insert(0,Z[6])
        self.Data.insert(0,Z)
    #----------------------------------------------- 
    def plt_R_I(self,Arg='o'):
        """Grafica la parte real e imaginária de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.subplot(211)     
        plt.semilogx(self.freqs[0],self.realZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$real\{Z\} [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.subplot(212)     
        plt.semilogx(self.freqs[0],self.imagZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$imag\{Z\} [\Omega]$')
        plt.show()
    def plt_Mag_Phase(self,Arg='o'):
        """Grafica la magnitud y fase de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.subplot(211)     
        plt.semilogx(self.freqs[0],self.magnZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$|Z| [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.subplot(212)     
        plt.semilogx(self.freqs[0],self.phasZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$Phase \{Z\} [Degrees]$')
        plt.show()
    def plt_Nyq(self,Arg='o'):
        """Grafica el diagrama de Nyquist"""
        import matplotlib.pyplot as plt
        plt.plot(self.realZ[0],self.imagZ[0],Arg)
        plt.xlabel(r'$real\{Z\} [\Omega]$')
        plt.ylabel(r'$imag\{Z\} [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.show()
    def plt_real(self,Arg='o'):
        """Grafica la parte real de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.semilogx(self.freqs[0],self.realZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$real\{Z\} [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.show()    
    def plt_imag(self,Arg='o'):
        """Grafica la parte imaginaria de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.semilogx(self.freqs[0],self.imagZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$imag\{Z\} [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.show()    
    def plt_magn(self,Arg='o'):
        """Grafica la magnitud de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.semilogx(self.freqs[0],self.magnZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$|Z| [\Omega]$')
        #plt.title('File name: ' + self.filename[0])
        plt.show()    
    def plt_phase(self,Arg='o'):
        """Grafica la fase de Z en función de 
        la frecuencia"""
        import matplotlib.pyplot as plt
        plt.semilogx(self.freqs[0],self.phasZ[0],Arg)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$phase\{Z\} [Degrees]$')
        #plt.title('File name: ' + self.filename[0])
        plt.show()
    def multiColefit(self,N=2,NIter=200,fMin=0,fMax=np.inf):
        """multiColefit ajusta los datos de impedancia al modelo de Cole-Cole
        de n dispersiones.  Los parámetros se almacenan como:
        x: array([Qi, DQ_1, tau_1, alpha_1, 
                        DQ_2, tau_2, alpha_2,
                        ... ,
                        DQ_n, tau_n, alpha_n])"""
        import scipy.optimize as opt
        import matplotlib.pyplot as plt
        f = self.freqs[0]
        Z = self.realZ[0] - 1j*self.imagZ[0]
        k = (f>=fMin) * (f<=fMax)
        f = f[k]
        Z = Z[k]
        w = 2 * np.pi * f
     
        kDQ    = range(1,3*N+1,3)
        ktau   = range(2,3*N+1,3)
        kalpha = range(3,3*N+1,3)
        
        def multicole(p,w1=w, k1=kDQ, k2=ktau, k3=kalpha):
            # multi Cole-Cole, Miranda et.al (2014)
            # p = [Qi,DQ_1,tau_1,alpha_1,DQ_2,tau_2,alpha_2, ..., DQ_n, tau_n,alpha_n]
            Qi    = p[0]
            DQ    = p[k1]
            tau   = p[k2]
            alpha = p[k3]
            M = w1.size
            N = DQ.size
            w1.resize(M,1)
            DQ.resize(1,N)
            tau.resize(1,N)
            alpha.resize(1,N)
            W    = w1 * np.ones((1,N))
            DQ   = np.ones((M,1))*DQ
            tau  = np.ones((M,1))*tau
            alpha= np.ones((M,1))*alpha
            aux = DQ/(1+(1j*W*tau)**(1-alpha))
            R = Qi + aux.sum(axis=1)
            return R
        
        wk = 1/Z.real - 1j/Z.imag
        
        p0    = np.zeros(1+3*N)
        p0[0] = (np.abs(Z.real.min())).tolist()
        aux2  = (1/w[(np.abs(Z.imag)).argmax()]).tolist()
        pDisp = [(Z.real.max()-Z.real.min())/N,aux2,1]
        for k in range(N):
            p0[1+3*k] = pDisp[0]
            p0[2+3*k] = pDisp[1]
            p0[3+3*k] = pDisp[2]

        def ffcole(p1,R0=Z,wk0=wk):
            if (p1 < 0).any():
                err = 100
            else:
                if p1[3] > 1:
                    err = 100
                else:
                    R1 = multicole(p1)
                    err = (((wk0.real*(R1.real - R0.real))**2 + 
                            (wk0.imag*(R1.imag - R0.imag))**2).sum())/R0.size
            return err
        #print 'Ajuste en proceso, por favor, espere ...'
        ps = opt.basinhopping(ffcole,p0,niter=NIter)
        
        f_fit = np.logspace(np.log10(self.freqs[0].min()),
                np.log10(self.freqs[0].max()),1000)
        w_fit = 2*np.pi*f_fit
        plt.subplot(211)
        h = plt.semilogx(f,Z.real,'o',
                    f_fit,multicole(ps.x,w_fit).real)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$real\{Z\} [\Omega]$')
        plt.legend(h,('Datos','Ajuste: Cole'))

        plt.subplot(212)     
        plt.semilogx(f,-Z.imag,'o',
                    f_fit,-multicole(ps.x,w_fit).imag)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$-imag\{Z\} [\Omega]$')
        plt.show()
        
        ks = ps.x[ktau].argsort()
        
        ps['Qi']      = ps.x[0]
        ps['Qo']      = ps.x[0] + ps.x[kDQ][ks].sum()
        ps['DQ_l']    = ps.x[kDQ][ks]                        
        ps['tau_l']   = ps.x[ktau][ks]
        ps['alpha_l'] = ps.x[kalpha][ks]
        self.p_multiColefit.insert(0,ps)
        return ps
    def p_multiColePlot(self,opt=0):
        """ Esta función grafica el ajuste para la optimización opt.
            Note que la última optimización realiada se almacena en
            opt=0"""
        Message = 'Usted ha realizado %d optimizaciones.\nSeleccione un valor apropiado para opt' % len(self.p_multiColefit)
        import matplotlib.pyplot as plt
        if opt+1 > len(self.p_multiColefit):
            p = Message
            print(Message)
        else:
            p = self.p_multiColefit[opt]
            f = self.freqs[0]
            Z = self.realZ[0] - 1j*self.imagZ[0]
            w  = 2 * np.pi * f
            kDQ    = range(1,len(p.x),3)
            ktau   = range(2,len(p.x),3)
            kalpha = range(3,len(p.x),3)
            def multicole(p,w1=w, k1=kDQ, k2=ktau, k3=kalpha):
                # multi Cole-Cole, Miranda et.al (2014)
                # p = [Qi,DQ_1,tau_1,alpha_1,DQ_2,tau_2,alpha_2, ..., DQ_n, tau_n,alpha_n]
                Qi    = p[0]
                DQ    = p[k1]
                tau   = p[k2]
                alpha = p[k3]
                M = w1.size
                N = DQ.size
                w1.resize(M,1)
                DQ.resize(1,N)
                tau.resize(1,N)
                alpha.resize(1,N)
                W    = w1 * np.ones((1,N))
                DQ   = np.ones((M,1))*DQ
                tau  = np.ones((M,1))*tau
                alpha= np.ones((M,1))*alpha
                aux = DQ/(1+(1j*W*tau)**(1-alpha))
                R = Qi + aux.sum(axis=1)
                return R
            f_fit = np.logspace(np.log10(self.freqs[0].min()),
                    np.log10(self.freqs[0].max()),1000)
            w_fit = 2*np.pi*f_fit
            plt.subplot(211)
            h = plt.semilogx(f,Z.real,'o',
                        f_fit,multicole(p.x,w_fit).real)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$real\{Z\} [\Omega]$')
            plt.legend(h,('Datos','Ajuste: MultiCole'))
    
            plt.subplot(212)     
            plt.semilogx(f,-Z.imag,'o',
                        f_fit,-multicole(p.x,w_fit).imag)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$-imag\{Z\} [\Omega]$')
            plt.show()
        return p
    def Colefit(self,NIter=200,fMin=0,fMax=np.inf):
        """Colefit ajusta los datos de impedancia al modelo de Cole-Cole
        de una dispersión.  Los parámetros se almacenan como:
        x: array([R0,Ri,tau,alpha])"""
        p = self.multiColefit(1,NIter,fMin,fMax)
        ps = p.copy()
        ps.pop('Qo');ps.pop('Qi');ps.pop('DQ_l');
        ps.pop('tau_l');ps.pop('alpha_l');ps.pop('x');
        ps['Ro']=p['Qo']
        ps['Ri']=p['Qi']
        ps['tau']=p['tau_l'][0]
        ps['alpha']=p['alpha_l'][0]
        ps['x']=[p['Qo'],p['Qi'],p['tau_l'][0],p['alpha_l'][0]]
        self.p_Colefit.insert(0,ps)
        return ps
    def p_ColePlot(self,opt=0):
        """ Esta función grafica el ajuste para la optimización opt.
            Note que la última optimización realiada se almacena en
            opt=0"""
        Message = 'Usted ha realizado %d optimizaciones.\nSeleccione un valor apropiado para opt' % len(self.p_multiColefit)
        import matplotlib.pyplot as plt
        if opt+1 > len(self.p_Colefit):
            p = Message
            print(Message)
        else:
            p = self.p_Colefit[opt]
            f = self.freqs[0]
            Z = self.realZ[0] - 1j*self.imagZ[0]
            f_fit = np.logspace(np.log10(self.freqs[0].min()),
                    np.log10(self.freqs[0].max()),1000)
            w_fit = 2*np.pi*f_fit
            Ro    = p['Ro']
            Ri    = p['Ri']
            tau   = p['tau']
            alpha = p['alpha']
            Z_fit = Ri + (Ro-Ri)/(1+(1j*w_fit*tau)**(1-alpha))
            plt.subplot(211)
            h = plt.semilogx(f,Z.real,'o',
                        f_fit,Z_fit.real)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$real\{Z\} [\Omega]$')
            plt.legend(h,('Datos','Ajuste: Cole-Cole'))
    
            plt.subplot(212)     
            plt.semilogx(f,-Z.imag,'o',
                        f_fit,-Z_fit.imag)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$-imag\{Z\} [\Omega]$')
            plt.show()
        return p
    def GEMTIP_sph(self,N=2,NIter=200,fMin=0,fMax=np.inf):
        """GEMTIP_sph ajusta los datos de impedancia al modelo
                de inclusiones esféricas tratadas con la teoría
                de medio efectivo, GEMTIP, para n dispersiones.
                Los parámetros se almacenan como:
        x: array([Ro, M:1, tau_1, c_1, 
                        M_2, tau_2, c_2,
                        ... ,
                        M_n, tau_n, c_n])"""
        import scipy.optimize as opt
        import matplotlib.pyplot as plt
        f = self.freqs[0]
        Z = self.realZ[0] - 1j*self.imagZ[0]
        Z0= Z.real.max()
        Z = Z / Z0
        k = (f>=fMin) * (f<=fMax)
        f = f[k]
        Z = Z[k]
        w = 2 * np.pi * f

        kMl  = range(1,3*N+1,3)
        ktau = range(2,3*N+1,3)
        kc   = range(3,3*N+1,3)
        
        def gemtip_sph(p,w1=w,k1=kMl, k2=ktau, k3=kc):
            # multi Cole-Cole, Miranda et.al (2014)
            # p = [Qi,DQ_1,tau_1,alpha_1,DQ_2,tau_2,alpha_2, ..., DQ_n, tau_n,alpha_n]
            Ro  = p[0]
            Ml  = p[k1]
            tau = p[k2]
            c   = p[k3]
            M   = w1.size
            N   = Ml.size
            w1.resize(M,1)
            Ml.resize(1,N)
            tau.resize(1,N)
            c.resize(1,N)
            W   = w1 * np.ones((1,N))
            Ml  = np.ones((M,1))*Ml
            tau = np.ones((M,1))*tau
            c   = np.ones((M,1))*c
            aux = Ml*( 1 - 1/( 1 + (1j*W*tau)**c ) )
            R = Ro/( 1 + aux.sum(axis=1))
            return R

        wk = 1/Z.real - 1j/Z.imag
        
        p0    = np.zeros(1+3*N)
        p0[0] = (Z.real.max()).tolist()
        aux2  = (1/w[(np.abs(Z.imag)).argmax()]).tolist()
        pDisp = [np.abs(Z.real.max()/Z.real.min() - 1)/N,aux2,0.5]
        for k in range(N):
            p0[1+3*k] = pDisp[0]
            p0[2+3*k] = pDisp[1]
            p0[3+3*k] = pDisp[2]
            
        def ffgsh(p1,R0=Z,wk0=wk,k1=kMl, k2=ktau, k3=kc):
            if (p1 < 0).any():
                err = 100
            else:
                if (p1[k1] > 3).any():
                    err = 100
                else:
                    if (p1[k3] > 1).any():
                        err = 100
                    else:
                        R1 = gemtip_sph(p1)
                        err = (((wk0.real*(R1.real - R0.real))**2 + 
                                (wk0.imag*(R1.imag - R0.imag))**2).sum())/R0.size
            return err
        
        print('Ajuste en proceso, por favor, espere ...')
        ps = opt.basinhopping(ffgsh,p0,niter=NIter)
        ps.x[0] = Z0 * ps.x[0]

        f_fit = np.logspace(np.log10(self.freqs[0].min()),
                np.log10(self.freqs[0].max()),1000)
        w_fit = 2*np.pi*f_fit
        plt.subplot(211)
        h = plt.semilogx(f,Z0*Z.real,'o',
                    f_fit,gemtip_sph(ps.x,w_fit).real)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$real\{Z\} [\Omega]$')
        plt.legend(h,('Datos','Ajuste: GEMTIP_sph'))

        plt.subplot(212)     
        plt.semilogx(f,-Z0*Z.imag,'o',
                    f_fit,-gemtip_sph(ps.x,w_fit).imag)
        plt.xlabel(r'$f [Hz]$')
        plt.ylabel(r'$-imag\{Z\} [\Omega]$')
        plt.show()
        ks = ps.x[ktau].argsort()
        ps['c_l']   = ps.x[kc][ks]
        ps['tau_l'] = ps.x[ktau][ks]
        ps['M_l']   = ps.x[kMl][ks]
        ps['Ri']    = ps.x[0] / (1 + ps.x[kMl][ks].sum())
        ps['Ro']    = ps.x[0]
        self.p_GEMTIPsphfit.insert(0,ps)
        return ps
    def p_GEMTIPsphPlot(self,opt=0):
        """ Esta función grafica el ajuste para la optimización opt.
            Note que la última optimización realiada se almacena en
            opt=0"""
        Message = 'Usted ha realizado %d optimizaciones.\nSeleccione un valor apropiado para opt' % len(self.p_multiColefit)
        import matplotlib.pyplot as plt
        if opt+1 > len(self.p_GEMTIPsphfit):
            p = Message
            print(Message)
        else:
            p = self.p_GEMTIPsphfit[opt]
            f = self.freqs[0]
            Z = self.realZ[0] - 1j*self.imagZ[0]
            f_fit = np.logspace(np.log10(self.freqs[0].min()),
                    np.log10(self.freqs[0].max()),1000)
            w_fit = 2*np.pi*f_fit
            Ro  = p['Ro']
            Ml  = p['M_l']
            tau = p['tau_l']
            c   = p['c_l']
            M   = w_fit.size
            N   = Ml.size
            w_fit.resize(M,1)
            Ml.resize(1,N)
            tau.resize(1,N)
            c.resize(1,N)
            W   = w_fit.copy() * np.ones((1,N))
            Ml  = np.ones((M,1))*Ml
            tau = np.ones((M,1))*tau
            c   = np.ones((M,1))*c
            aux = Ml*( 1 - 1/( 1 + (1j*W*tau)**c ) )
            Z_fit = Ro/( 1 + aux.sum(axis=1))

            plt.subplot(211)
            h = plt.semilogx(f,Z.real,'o',
                        f_fit,Z_fit.real)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$real\{Z\} [\Omega]$')
            plt.legend(h,('Datos','Ajuste: GEMTIP_sph'))
            plt.subplot(212)     
            plt.semilogx(f,-Z.imag,'o',
                        f_fit,-Z_fit.imag)
            plt.xlabel(r'$f [Hz]$')
            plt.ylabel(r'$-imag\{Z\} [\Omega]$')
            plt.show()
        return p
def getzdata(dirname,filename,sep=';'):
    """Con esta función podrá obtener datos del programa AutoLab
    y genera una clase de Datos donde están los datos y funciones
    para graficar"""
    ZData = Datos(dirname=[], filename=[],
    Units=['[Z in Ohm], [frequencies in Hz], [Phase in degrees]'],
    Index=[],freqs=[],realZ=[],imagZ=[],magnZ=[],phasZ=[],Times=[],Data=[])
    #--- Lectura de datos ---
    f=open(os.path.join(dirname,filename),'r')
    aa = []; ka = 0
    for line in f:
        aa.insert(ka,line[:-1].split(sep))
        ka = ka + 1          
    f.close()
    # ------------------------
    Z = []
    Z = (np.matrix(aa[1:]).transpose()).tolist()
    Z = np.array(Z)
    
    ZData.rellenar(dirname,filename,aa[0],Z)  # Rellena los datos para cada paciente
    ZData.plt_R_I('*')
    return ZData

# getzdata('G:/Mi unidad/Colab_DanielTriana/data/20200622', 'Z_p56_1.txt')