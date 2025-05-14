from dustpy.simulation import Simulation
import numpy as np
import dustpy.constants as c
from wind import *
from simple_cal import *
from funcs import *
from fried_script import mdot_r, r_thin, sigma_dot_pe


sim = Simulation()

sim.ini.grid.Nr = 300
sim.ini.grid.rmin = 5* c.au
sim.ini.grid.rmax = 500 * c.au
sim.ini.star.M = 1 * c.M_sun
sim.ini.gas.SigmaRc =30 * c.au
sim.ini.gas.gamma = 1
sim.ini.gas.Mdisk= 0.02 *c.M_sun
sim.ini.gas.alpha_dw = 0
sim.ini.gas.alpha = 1e-3
sim.ini.dust.vfrag = 100
sim.ini.dust.d2gRatio= 0.01 # value to achieve Md/Mg=0.01

sim.ini.binary.M2= 0

sim.initialize()

sim.addgroup('fried', description='External PE from FRIED') #
sim.fried.addfield('pah', 0.1, description='PAH-to-dust mass ratio.')
sim.fried.addfield('dust', True, description='whether dust growth is considered.')
sim.fried.addfield('fuv', 1e3, description='strengths of FUV field (units: G0).')
sim.fried.addfield('mdot_pe', np.zeros(sim.grid.Nr), description='mass-loss ratio from FRIED (units: Msun/yr).')
sim.fried.addfield('rt_ind', 1, description='index of radius transiting to optically thin.' )

sim.fried.updater=['pah', 'dust', 'fuv', 'mdot_pe', 'rt_ind']
sim.updater=['star', 'grid', 'gas', 'dust', 'binary', 'fried']
sim.fried.mdot_pe.updater  = mdot_r
sim.fried.rt_ind.updater  = r_thin
sim.gas.S.ext.updater = sigma_dot_pe

gas_only(sim)
simple_save(sim)

sim.update()

snapshots = np.linspace(1e5, 3e6 , 30) * c.year
sim.writer.datadir = ''
sim.writer.overwrite = True
sim.t.snapshots = snapshots

sim.run()