"""
Quick Equilibration of Trypsin + 13U system (for testing)
Stage 1: NVT (heating 0 -> 300 K, 10 ps)
Stage 2: NPT (pressure, 10 ps)
Total: 20 ps (fast version)
"""

from openmm.app import *
from openmm import *
from openmm.unit import *

print("=" * 60)
print("EQUILIBRATION (FAST) - Trypsin + 13U + Water + Ions")
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

# ═══════════════════════════════════════════════════════════
# 2. SETUP PLATFORM (CPU)
# ═══════════════════════════════════════════════════════════
print("\nSTEP 2: Setting up platform...")

platform = Platform.getPlatformByName('CPU')
print(f"[OK] Platform: {platform.getName()}")

# ═══════════════════════════════════════════════════════════
# 3. LOAD MINIMIZED STATE
# ═══════════════════════════════════════════════════════════
print("\nSTEP 3: Loading minimized positions...")

with open('md/preparation/minimized_state.xml') as f:
    min_state = XmlSerializer.deserialize(f.read())

min_positions = min_state.getPositions()
print(f"[OK] Minimized positions loaded")

# ═══════════════════════════════════════════════════════════
# 4. STAGE 1: NVT (Heating 0 -> 300 K)
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("STAGE 1: NVT Equilibration (Heating)")
print("=" * 60)
print("  Temperature: 0 -> 300 K")
print("  Duration: 10 ps")
print("  Time (CPU): ~5-10 minutes")
print("=" * 60)

# Integrator with 2 fs (because HBonds constraints are in system)
integrator_nvt = LangevinMiddleIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds    # 2 fs (safe with HBonds)
)

simulation_nvt = Simulation(pdb.topology, system, integrator_nvt, platform)
simulation_nvt.context.setPositions(min_positions)

# Set initial velocities at 300 K
simulation_nvt.context.setVelocitiesToTemperature(300*kelvin)

print("[OK] NVT simulation created")
print("     Running 5,000 steps (10 ps)...")

# Report energy during equilibration
simulation_nvt.reporters.append(
    StateDataReporter(
        'md/preparation/nvt_log.csv',
        reportInterval=1000,
        step=True,
        temperature=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        speed=True
    )
)

# Run NVT
simulation_nvt.step(5000)  # 10 ps
print("[OK] Stage 1 (NVT) done")

# Save NVT state
state_nvt = simulation_nvt.context.getState(getPositions=True, getVelocities=True)
with open('md/preparation/nvt_state.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state_nvt))
print("[OK] nvt_state.xml saved")

# ═══════════════════════════════════════════════════════════
# 5. STAGE 2: NPT (Pressure 1 atm)
# ═══════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("STAGE 2: NPT Equilibration (Pressure)")
print("=" * 60)
print("  Temperature: 300 K")
print("  Pressure: 1 atm")
print("  Duration: 10 ps")
print("  Time (CPU): ~5-10 minutes")
print("=" * 60)

# Add barostat for NPT
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# New integrator for NPT
integrator_npt = LangevinMiddleIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds
)

simulation_npt = Simulation(pdb.topology, system, integrator_npt, platform)
simulation_npt.context.setPositions(state_nvt.getPositions())
simulation_npt.context.setVelocities(state_nvt.getVelocities())

# Report energy during NPT
simulation_npt.reporters.append(
    StateDataReporter(
        'md/preparation/npt_log.csv',
        reportInterval=1000,
        step=True,
        temperature=True,
        potentialEnergy=True,
        volume=True,
        density=True,
        speed=True
    )
)

print("[OK] NPT simulation created")
print("     Running 5,000 steps (10 ps)...")

# Run NPT
simulation_npt.step(5000)  # 10 ps
print("[OK] Stage 2 (NPT) done")

# Save final equilibrated state
state_final = simulation_npt.context.getState(getPositions=True, getVelocities=True)

with open('md/preparation/equilibrated_state.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state_final))
print("[OK] equilibrated_state.xml saved")

with open('md/preparation/equilibrated.pdb', 'w') as f:
    PDBFile.writeFile(pdb.topology, state_final.getPositions(), f)
print("[OK] equilibrated.pdb saved")

print("\n" + "=" * 60)
print("EQUILIBRATION COMPLETE!")
print("=" * 60)
print("  Total time: 20 ps")
print("  Final state ready for production MD")
print("=" * 60)