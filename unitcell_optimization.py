from ase.io import read
from ase.io import write
from ase.build import surface
from ase.visualize import view
from ase.constraints import *
from gpaw import GPAW
from ase.optimize import BFGS
from ase.io import Trajectory
from ase.constraints import StrainFilter
from gpaw import PW

#Include the unitcell structure
graphene = read("###.cif")

#Setting up GPAW as the calculator
#Select XC functional, I think the default is LDA
#We want to use Plane Waves with a wavefunction cutoff which 
# must be converged.
calc = GPAW(xc="PBE", mode=PW(400),
            kpts=(4,4,4))
graphene.set_calculator(calc)

#There are differnet methods to do unitcell optimization so check out the 
# GPAW website

sf = StrainFilter(graphene)

#Also check out other optimization algorithms on ASE!

#Each BFGS step with the calculated force will be written to the log file

opt = BFGS(sf, logfile="graphene.log")

#The optimized structure is saved in the .traj file
traj = Trajectory("graphene.traj", "w", graphene)
opt.attach(traj)

# fmax is the total force acting on each atom.
opt.run(fmax=0.025)
