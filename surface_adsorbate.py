from ase import Atoms
from ase.build import surface
from ase.build import add_adsorbate
from ase.io import read

# Build a slab from optimized structure
# For graphene it would the traj file

slab = read("###")

#Since you are dealing with grapehen you just expand it in x and y directions
# The procedure will be different for metal slabs with different number of layers and facets
slab = slab * (4,4,1)

#use Atoms class to build your adsorbate
# I recommend checking the Molecule class as well in ASE.
adsorbate = Atoms("2N", positions=[[0,0,0], [0,0,1.1]])

slab.set_calculator(calc)
E_slab = slab.get_potential_energy()
#Add adsorbate on the surface
# Check out the attributes of add_adsorbate function
#The default position is (0,0)
add_adsorbate(slab, adsorbate, 2)

calc = GPAW(mode=PW(400),
            xc='BEEF-vdW',
            kpts=(4, 4, 1))

# Perform the DFT calculation for the clean surface:
slab.set_calculator(calc)
E_total = slab.get_potential_energy()

# Perform the DFT calculation for the molecule adsorbed on the surface:
adsorbate.set_calculator(calc)
E_ads = adsorbate.get_potential_energy()

# Calculate the adsorption energy:
E_adsorption = E_total - E_slab - E_adsorbate



