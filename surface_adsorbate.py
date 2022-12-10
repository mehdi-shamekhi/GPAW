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

#Add adsorbate on the surface
# Check out the attributes of add_adsorbate function
#The default position is (0,0)
add_adsorbate(slab, adsorbate, 2)