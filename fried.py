import numpy as np
import pandas as pd
import math
from scipy.interpolate import LinearNDInterpolator, interp1d
import os

class FRIED():
        
    def __init__(self, mstar=0.1, pah=0.1, dust=True, fuv=10, Nr= 1000, Ns=1000, radial=None, sigma_g=None):

        mstar_ls = [0.1, 0.3, 0.6, 1.0, 1.5, 3.0]
        pah_ls   = [0.1, 0.5]
        fuv_ls   = [1, 10, 1e2, 1e3, 1e4, 1e5]

        if mstar not in mstar_ls:
            raise Exception("Stellar mass should be selected from 0.1, 0.3, 0.6, 1.0, 1.5, 3.0")
        if pah not in pah_ls:
            raise Exception("The value of PAH should be selected from 0.1, 0.5")
        if fuv not in fuv_ls:
            raise Exception("The FUV stength should be selected from 1, 10, 100, 1000, 1e4, 1e5")


        self.mstar = mstar
        self.pah_ratio = pah
        self.dust_grow = dust
        self.FUV       = fuv

        # interpolated grid
        self.Nr        = Nr # radial grid
        self.Ns        = Ns # gas surface density

        if self.dust_grow == True:
            suffix = 'growth'
        else:
            suffix = 'ism'

        frac_m, whole_m = math.modf(self.mstar)
        frac_m = int(frac_m*10)
        whole_m = int(whole_m)
        self.mstar = f'{whole_m}p{frac_m}'

        frac_p, whole_p = math.modf(self.pah_ratio)
        frac_p = int(frac_p*10)
        whole_p = int(whole_p)
        self.pah_ratio  = f'{whole_p}p{frac_p}'
        
        self.name = f'FRIEDV2_{self.mstar}Msol_fPAH{self.pah_ratio}_{suffix}.dat'
        self.data      = self.read_data()
        self.sigma_uni, self.N, self.M = self.data_shape()
        self.sigma_hydro = np.zeros((self.N,self.M))

        self.r_new     = None
        self.sigma_new = None
        self.mdot_new  = None
        
        self.radial    = radial
        self.sigma_g   = sigma_g
        
        
    def read_data(self):
        
        path=os.getcwd()
        df = pd.read_csv(path+'/data/'+self.name, names = ['Host star mass [Msol]', 'Disc outer radius [au]', 
                                'Surface density at 1au [gcm-2]', 'Surface density at disc outer edge [gcm-2]',
                                'FUV field strength at outer edge of grid [G0]', 'Mass loss rate [log10(Msol/yr)]'])
        df2 = df[df['FUV field strength at outer edge of grid [G0]'] == self.FUV]

        #self.data = df2
        
        return df2

    def data_shape(self):

        '''
        compute how many gas surface densities do we use
        in hydrodynamic simulations (N); and how many samples 
        do we have for each gas surface density (M). 
         '''

        self.sigma_uni = np.unique(self.data['Surface density at 1au [gcm-2]'].values)
        self.N =  np.unique(self.data['Surface density at 1au [gcm-2]'].values).shape[0]# how many surface densities do we have
        self.M = (self.data[self.data['Surface density at 1au [gcm-2]']==self.sigma_uni[0]]).shape[0] # samples we have for one specific surface density

        return self.sigma_uni, self.N, self.M

    def hydro_sigma_g(self):

        '''collect gas surface densities used in hydro simulations

        (multiple sampling for one specific surface density)
        
        return: surface density in an array (self.N*self.M)
        
        '''

        for ind, i in enumerate(self.sigma_uni):
            data2 = self.data[self.data['Surface density at 1au [gcm-2]'].values==i]
            self.sigma_hydro[ind, :]= data2['Surface density at disc outer edge [gcm-2]'].values

        return self.sigma_hydro

    def interp_hdyro_mdot(self):

        '''collect the mass-loss rate from hydro simulations
        and interpolate it to the new grid.
        '''
        
        
        r = self.data['Disc outer radius [au]'].values
        sigma = self.data['Surface density at disc outer edge [gcm-2]'].values
        mdot = 10**(self.data['Mass loss rate [log10(Msol/yr)]'].values)
        interp = LinearNDInterpolator(list(zip(r, sigma)), mdot)

        # interpolate to new grid
        if (self.radial is None) | (self.sigma_g is None):
            print('No radial grid/ gas surface density is given.')
            print('Using the default grid.') 
            self.r_new = np.linspace(5, r.max(), self.Nr)[::-1]
            self.sigma_new = np.logspace(np.log10(sigma.min()), np.log10(sigma.max()), self.Ns)
        #    
        else:
            self.r_new = self.radial
            self.sigma_new = self.sigma_g
            #xx  = self.r_new
            #yy  = self.sigma_new
        
        xx, yy = np.meshgrid(self.r_new, self.sigma_new)
        mdot_new = interp(xx, yy).T
       
        # remove mdot estimated from intepolation beyond the hydro data.
        r_uni = np.unique(r)
        sigma_ = self.hydro_sigma_g()

        # interpolate the gas surface density in hydro sim 
        # to find the upper & lower limits of gas surface density 
        # at each new radial cell

        sigma_up_lim_fit = interp1d(r_uni, sigma_[-1,:][::-1])
        sigma_lo_lim_fit = interp1d(r_uni, sigma_[0,:][::-1])

        sigma_up = sigma_up_lim_fit(self.r_new)
        sigma_lo = sigma_lo_lim_fit(self.r_new)

        for ind, val in enumerate(self.r_new):
            for ind_2, val2 in enumerate(self.sigma_new):
                if (val2> sigma_up[ind]) | (val2 < sigma_lo[ind]):
                    mdot_new[ind, ind_2] = np.nan
        
        self.mdot_new = mdot_new

        return

    def extrap_hydro_mdot(self):

        for ind, val in enumerate(self.r_new):

            loc_not_nan = np.where(~np.isnan(self.mdot_new[ind,:]))
            loc_nan = np.where(np.isnan(self.mdot_new[ind,:])) # location with (no) interpolate mdot
            if len(loc_not_nan[0])!=0:
                sigma_max = self.sigma_new[loc_not_nan].max()
                sigma_min = self.sigma_new[loc_not_nan].min()
                mdot_max  = np.nanmax(self.mdot_new[ind, :])
                mdot_min  = np.nanmin(self.mdot_new[ind, :])
        
            for i in loc_nan[0]:
                if self.sigma_new[i]>sigma_max:
                    self.mdot_new[ind, i]=mdot_max
                elif self.sigma_new[i]<sigma_min:
                    self.mdot_new[ind, i]= mdot_min * self.sigma_new[i]/sigma_min
        
            else:
                pass

        return

    
    #def mdot_r(self):
    #    '''compute the mass-loss rate as a function of radius from FRIED grid
    #    
    #    Input:
    #    ------
    #    sim from DustPy for radial grid and gas surface density
    #    mstar: mass of central star [Msun]
    #    pah: ratio of PAH to dust ??
    #    dust: boolean, whether consider dust growth
    #    fuv : strength of FUV
    #
    #    Output:
    #    -------
    #    mdot: Mass-loss rate as a function of radius [Msun/year]
    #    '''
    #
    #    f = FRIED(self.mstar, self.pah_ratio, self.dust_grow, self.FUV, self.radial, self.sigma_g)
    #    f.interp_hdyro_mdot()
    #    f.extrap_hydro_mdot()
    #    mdot= np.zeros(len(sim.grid.r))
    #    for i in range(len(sim.grid.r)):
    #        mdot[i] = f.mdot_new[i, i]
    #
    #    return mdot


        