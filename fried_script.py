from FRIEDgrid.fried import FRIED # server/Jupyter Notebook
#from fried import FRIED # server/jobs
import dustpy.constants as c
from scipy.interpolate import interp1d
import numpy as np
import astropy.constants as astro_c
from astropy import units as u

# read the 

def mdot_r(sim):
    '''
    Return the interpolated and extrapolated mass-loss rate as a
    function of radius for the given radial grid.

    Input: 
    ------
    mdot: mass loss rate as a function of radius

    '''
    f = FRIED(sim.star.M/c.M_sun, sim.fried.pah, sim.fried.dust, sim.fried.fuv, radial=sim.grid.r/c.au, sigma_g=sim.gas.Sigma)
    f.interp_hdyro_mdot()
    f.extrap_hydro_mdot()
    
    mdot = np.zeros(len(sim.grid.r))
    for i in range(sim.grid.Nr):
        mdot[i] = f.mdot_new[i,i] 

    # convert Msun/yr to cgs
    mdot = mdot * astro_c.M_sun.cgs.value /u.year.to(u.second)
    return mdot

def r_thin(sim):

    '''identify the radius where the gas disc becomes optically thin;
    
    Input: 
    ------
    sim: from DustPy  
    mdot: the mass-loss rate read from FRIED grid.

    Output:
    -------
    the index of the radius transiting from optically thick to thin discs
    the index of the transition radius
    '''
    mdot= sim.fried.mdot_pe
    # second derivative of the mass-loss rate
    face =interp1d(sim.grid.r, mdot, fill_value='extrapolate') # mdot at grid face
    mdot_face = face(sim.grid.ri)
    grad = np.diff(mdot_face)/np.diff(sim.grid.ri)  
    face2 = interp1d(sim.grid.r, grad, fill_value='extrapolate')
    grad_face = face2(sim.grid.ri)
    grad2nd = np.abs(np.diff(grad_face)/np.diff(sim.grid.ri))

    rt_ind= np.argmax(grad2nd)
    #r_t  = sim.grid.r[rt_ind]

    return rt_ind

def sigma_dot_pe(sim):

    '''Compute the effective loss rate of each grid in the 
    optically thin disc [in units of surface density]; this
    can be directly applied to DustPy as the external mass loss
    
    Eq. 6 & 7 in Sellek et al. (2020)
    '''

    rt_ind = sim.fried.rt_ind
    mdot   = sim.fried.mdot_pe

    #Eq. 6 in Sellek et al. (2020)

    dr = sim.grid.ri[rt_ind+1:]-sim.grid.ri[rt_ind:-1]
    r  = sim.grid.r[rt_ind:]
    sigma_g_ = sim.gas.Sigma[rt_ind:]

    Mi = 2*np.pi*r*dr*sigma_g_
    M_beyond = np.sum(Mi)

    mdot_tot = np.sum(Mi/M_beyond * mdot[rt_ind:])

    # Eq. 7 in Sellek et al. (2020)
    sigma_rt  = np.copy(sim.gas.Sigma)
    np.put(sigma_rt, np.arange(rt_ind), np.zeros(rt_ind)) 
    sigma_dot = -1*sigma_rt * mdot_tot/M_beyond

    return sigma_dot