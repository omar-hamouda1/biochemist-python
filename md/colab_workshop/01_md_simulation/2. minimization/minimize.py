"""
Minimization of Trypsin + 13U system
Goal: Remove clashes before MD
Two-stage approach with restraints (safe settings)
"""

from openmm.app import *
from openmm import *
from openmm.unit import *

print("=" * 60)
print("MINIMIZATION - Trypsin + 13U + Water + Ions")
print("=" * 60)

# ═══════════════════════════════════════════════════════════
# 1. LOAD SYSTEM
# ═══════════════════════════════════════════════════════════
print("\nSTEP 1: Loading system...")

pdb = PDBFile('md/preparation/system.pdb')

with open('md/preparation/system.xml') as f:
    system = XmlSerializer.deserialize(f.read())

print(f"[OK] System loaded")
print(f"     Atoms: {system.getNumParticles()}")
print(f"     Forces: {system.getNumForces()}")

# ═══════════════════════════════════════════════════════════
# 2. SETUP PLATFORM (CPU)
# ═══════════════════════════════════════════════════════════
print("\nSTEP 2: Setting up platform...")

platform = Platform.getPlatformByName('CPU')
print(f"[OK] Platform: {platform.getName()}")

# ═══════════════════════════════════════════════════════════
# 3. STAGE 1: MINIMIZE WITH RESTRAINTS ON HEAVY ATOMS
# ═══════════════════════════════════════════════════════════
print("\nSTEP 3: Stage 1 - Minimize with restraints on heavy atoms...")

# Add position restraints on protein heavy atoms (LIGHT restraints)
restraint = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
restraint.addGlobalParameter("k", 100.0)   # REDUCED from 1000 to 100
restraint.addPerParticleParameter("x0")
restraint.addPerParticleParameter("y0")
restraint.addPerParticleParameter("z0")

positions = pdb.positions
restrained_count = 0

for atom in pdb.topology.atoms():
    if atom.element.symbol != 'H' and atom.residue.name not in ['HOH', 'WAT', 'NA', 'CL']:
        pos = positions[atom.index]
        restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
        restrained_count += 1

system.addForce(restraint)
print(f"[OK] Restrained {restrained_count} heavy atoms (k=100)")

# NEW integrator for stage 1
integrator1 = LangevinMiddleIntegrator(
    300*kelvin,
    1/picosecond,
    0.001*picoseconds
)

simulation = Simulation(pdb.topology, system, integrator1, platform)
simulation.context.setPositions(pdb.positions)

# Calculate initial energy
state = simulation.context.getState(getEnergy=True)
initial_energy = state.getPotentialEnergy()
print(f"[!] Initial energy: {initial_energy}")

# Minimize stage 1 (SAFE settings)
print("     Running stage 1 (with restraints)...")
simulation.minimizeEnergy(
    tolerance=1000*kilojoules_per_mole/nanometer,   # LARGER tolerance
    maxIterations=1000                               # FEWER iterations
)
print("[OK] Stage 1 done")

# Get positions after stage 1
state = simulation.context.getState(getPositions=True)
pos_stage1 = state.getPositions()

# ═══════════════════════════════════════════════════════════
# 4. STAGE 2: MINIMIZE WITHOUT RESTRAINTS
# ═══════════════════════════════════════════════════════════
print("\nSTEP 4: Stage 2 - Minimize without restraints...")

# Remove restraint force
system.removeForce(system.getNumForces() - 1)

# NEW integrator for stage 2
integrator2 = LangevinMiddleIntegrator(
    300*kelvin,
    1/picosecond,
    0.001*picoseconds
)

# New simulation without restraints
simulation2 = Simulation(pdb.topology, system, integrator2, platform)
simulation2.context.setPositions(pos_stage1)

# Minimize stage 2
print("     Running stage 2 (no restraints)...")
simulation2.minimizeEnergy(
    tolerance=100*kilojoules_per_mole/nanometer,
    maxIterations=2000
)
print("[OK] Stage 2 done")

# ═══════════════════════════════════════════════════════════
# 5. FINAL ENERGY
# ═══════════════════════════════════════════════════════════
print("\nSTEP 5: Calculating final energy...")

state = simulation2.context.getState(getEnergy=True, getPositions=True)
final_energy = state.getPotentialEnergy()
print(f"[OK] Final energy: {final_energy}")

improvement = (initial_energy - final_energy) / initial_energy * 100
print(f"\n[STATS] IMPROVEMENT: {improvement:.2f}%")

# ═══════════════════════════════════════════════════════════
# 6. SAVE
# ═══════════════════════════════════════════════════════════
print("\nSTEP 6: Saving minimized structure...")

positions = state.getPositions()
with open('md/preparation/minimized.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, positions, f)
print("[OK] minimized.pdb saved")

with open('md/preparation/minimized_state.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
print("[OK] minimized_state.xml saved")

print("\n" + "=" * 60)
print("MINIMIZATION COMPLETE!")
print("=" * 60)