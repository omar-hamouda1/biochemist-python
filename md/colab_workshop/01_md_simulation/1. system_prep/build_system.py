"""
Build MD System: Trypsin + 13U ligand + water + ions
Using RDKit -> OpenFF -> GAFFTemplateGenerator
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit import Molecule

print("=" * 60)
print("SYSTEM BUILDING - Trypsin + 13U + Water + Ions")
print("=" * 60)

# ═══════════════════════════════════════════════════════════
# 1. LOAD PROTEIN
# ═══════════════════════════════════════════════════════════
print("\nSTEP 1: Loading protein...")
protein_pdb = PDBFile('pdb/protein_h.pdb')
print(f"[OK] Protein loaded: {protein_pdb.topology.getNumAtoms()} atoms")

# ═══════════════════════════════════════════════════════════
# 2. LOAD LIGAND WITH RDKIT
# ═══════════════════════════════════════════════════════════
print("\nSTEP 2: Loading ligand with RDKit...")
rdkit_mol = Chem.MolFromMolFile('docking/results/fixed/13U_fixed.sdf')
rdkit_mol = Chem.AddHs(rdkit_mol)
print(f"[OK] Ligand: {rdkit_mol.GetNumAtoms()} atoms (with H)")

# ═══════════════════════════════════════════════════════════
# 3. CONVERT RDKit -> OpenFF
# ═══════════════════════════════════════════════════════════
print("\nSTEP 3: Converting RDKit -> OpenFF...")
ligand_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
print(f"[OK] OpenFF Molecule: {ligand_mol.n_atoms} atoms")

# ═══════════════════════════════════════════════════════════
# 4. CREATE GAFF TEMPLATE GENERATOR (AM1-BCC)
# ═══════════════════════════════════════════════════════════
print("\nSTEP 4: Creating GAFFTemplateGenerator (AM1-BCC)...")
print("     This will compute AM1-BCC charges. May take 1-3 minutes...")

gaff = GAFFTemplateGenerator(
    molecules=ligand_mol,
    forcefield='gaff-2.11'
)
print("[OK] GAFFTemplateGenerator created")

# ═══════════════════════════════════════════════════════════
# 5. CREATE FORCE FIELD WITH GAFF
# ═══════════════════════════════════════════════════════════
print("\nSTEP 5: Creating ForceField...")
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
forcefield.registerTemplateGenerator(gaff.generator)
print("[OK] ForceField created with GAFF")

# ═══════════════════════════════════════════════════════════
# 6. CREATE MODELLER AND ADD LIGAND
# ═══════════════════════════════════════════════════════════
print("\nSTEP 6: Creating Modeller and adding ligand...")

ligand_topology = ligand_mol.to_topology().to_openmm()
ligand_positions = ligand_mol.conformers[0].to_openmm()

modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.add(ligand_topology, ligand_positions)
print(f"[OK] Complex: {modeller.topology.getNumAtoms()} atoms")

# ═══════════════════════════════════════════════════════════
# 7. ADD WATER + IONS
# ═══════════════════════════════════════════════════════════
print("\nSTEP 7: Adding water + ions...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0*nanometer,
    ionicStrength=0.15*molar,
    neutralize=True
)
print(f"[OK] Solvent added: {modeller.topology.getNumAtoms()} atoms")

# ═══════════════════════════════════════════════════════════
# 8. CREATE SYSTEM
# ═══════════════════════════════════════════════════════════
print("\nSTEP 8: Creating OpenMM System...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.0*nanometer,
    constraints=HBonds,
    rigidWater=True
)
print(f"[OK] System created: {system.getNumParticles()} particles")

# ═══════════════════════════════════════════════════════════
# 9. SAVE SYSTEM
# ═══════════════════════════════════════════════════════════
print("\nSTEP 9: Saving system files...")

with open('md/preparation/system.pdb', 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("[OK] system.pdb saved")

with open('md/preparation/system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(system))
print("[OK] system.xml saved")

print("\n" + "=" * 60)
print("SYSTEM BUILDING COMPLETE!")
print("=" * 60)