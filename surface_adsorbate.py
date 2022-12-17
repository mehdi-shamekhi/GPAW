from ase import Atoms
from ase.build import surface
from ase.build import add_adsorbate
from ase.io import read
from gpaw import GPAW, PW
from ase.optimize import BFGS

#Read the optimized unitcell

unitcell = read("./C.cif")

slab = surface(unitcell, (0,0,1), 1, vacuum=10)

calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(1, 1, 1))

slab.set_calculator(calc)

opt = BFGS(slab, trajectory='slab.traj',
           logfile='slab.log')

opt.run(fmax=0.05)

E_slab = slab.get_potential_energy()

adsorbate = Atoms("2N", positions=[[0,0,0], [0,0,1.1]])

add_adsorbate(slab, adsorbate, 2)

opt = BFGS(slab, trajectory='slab_adsorbate.traj',
           logfile='slab_adsorbate.log')
opt.run(fmax=0.05)

E_total = slab.get_potential_energy()


adsorbate = Atoms("2N", positions=[[0,0,0], [0,0,1.1]], cell=[10,10,10])
adsorbate.center()

adsorbate.set_calculator(calc)

opt = BFGS(adsorbate, trajectory='adsorbate.traj',
           logfile='adsorbate.log')
opt.run(fmax=0.05)


E_adsorption = E_total - E_slab - E_adsorbate

print(E_adsorption)
