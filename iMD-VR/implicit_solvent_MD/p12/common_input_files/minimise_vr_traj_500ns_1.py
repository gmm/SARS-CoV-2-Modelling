from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md
import numpy as np
import math
from sys import stdout


#4 minimization steps: 5 kcal/mol restraints, half after each. Run at 10K. 
#At around 1 kcal/mol, only restrain bbone alpha carbons and ligand. No restraints on last step.

#heat over 200ps to 298K, restrain backbone

#8 equilibration steps, of around 500ps? each. 5 kcal/mol 
#restraint, half each time. reduce to bbone alpha carbons and ligand. No restraints on last step.



def make_new_system(prmtop, restraint_atoms, restraint_force, positions):
    system = prmtop.createSystem(nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=0.8 * nanometer, constraints=None,
                                 implicitSolvent=OBC2)

    force = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    force.addGlobalParameter('k', restraint_force * kilocalories_per_mole / angstroms ** 2)
    force.addPerParticleParameter('x0')
    force.addPerParticleParameter('y0')
    force.addPerParticleParameter('z0')
    for i in restraint_atoms:
        force.addParticle(int(i), positions[int(i)].value_in_unit(nanometers))
    system.addForce(force)

    return system


pdb = 'vr_docking/min_hbonds_1.pdb'
prmtop = AmberPrmtopFile('../make_apo_subs_prmtop/6yb7_and_pep12.prmtop')

### Define restraints #####

print('the name of the input pdb is:', pdb)

# All atoms
all_atoms =  md.Trajectory.load(pdb).topology.select('protein')

# Main protease atoms
mpro_atoms = md.Trajectory.load(pdb).topology.select('chainid 0 1')

# Backbone atoms
bbone_atoms = md.Trajectory.load(pdb).topology.select('backbone and chainid 0 1')

# Polypeptide substrate
subs_atoms = md.Trajectory.load(pdb).topology.select('chainid 2')

# Alpha atoms
alpha_atoms = md.Trajectory.load(pdb).topology.select('chainid 0 1 and name CA')

print('the total number of atoms is:', len(all_atoms))
print('the number of atoms in the main protease is:', len(mpro_atoms))
print('the number of substrate atoms is:', len(subs_atoms))
print('the number of atoms in the backbone is:', len(bbone_atoms))
print('the number of alpha carbon atoms is:', len(alpha_atoms))

### Minimization #####

# Set up simulation
pdb = PDBFile(pdb)
positions = pdb.positions
platform = Platform.getPlatformByName('OpenCL')

# Minimize
rounds = 4
restraint_force = 5
restraint_atoms = bbone_atoms

for i in range(0, rounds):
    system = make_new_system(prmtop, restraint_atoms, restraint_force, positions)
    integrator = LangevinIntegrator(10 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    print('Minimizing: Restraint of ' + str(restraint_force) + ' kcal/mol, applied to ' + str(
        len(restraint_atoms)) + ' atoms.')
    simulation.minimizeEnergy()

    if i == (rounds - 2):
        restraint_atoms = []
        restraint_force = 0
    else:
        restraint_force = restraint_force / 2
        if restraint_force < 1:
            restraint_atoms = alpha_atoms

    positions = simulation.context.getState(getPositions=True).getPositions()

### Heating #####

restraint_atoms = bbone_atoms
system = make_new_system(prmtop, restraint_atoms, 5, positions)

heating_increments = [10, 40, 70, 100, 130, 160, 190, 220, 250, 298]

total_sim_time = 0

for inc in heating_increments:
    print('Heating to ' + str(inc) + ' K...')
    integrator = LangevinIntegrator(inc * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
                                                      potentialEnergy=True, temperature=True, progress=True,
                                                      remainingTime=True,
                                                      speed=True, totalSteps=1000, separator='\t'))
    simulation.reporters.append(DCDReporter(f'binding_1_production_md_500ns_1/heat_{inc}_prod_6yb7_pep12.dcd', 1000))
    simulation.step(10000)
    total_sim_time += 10000 * 0.002

    positions = simulation.context.getState(getPositions=True).getPositions()

print(total_sim_time)

### Equilibration #####
rounds = 8
restraint_force = 5
restraint_atoms = bbone_atoms

for i in range(0, rounds):
    system = make_new_system(prmtop, restraint_atoms, restraint_force, positions)
    integrator = LangevinIntegrator(298 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(positions)
    print('Equilibrating: Restraint of ' + str(restraint_force) + ' kcal/mol, applied to ' + str(
        len(restraint_atoms)) + ' atoms.')
    simulation.reporters.append(app.StateDataReporter(stdout, 5000, step=True,
                                                      potentialEnergy=True, temperature=True, progress=True,
                                                      remainingTime=True,
                                                      speed=True, totalSteps=1000, separator='\t'))
    simulation.reporters.append(DCDReporter(f'binding_1_production_md_500ns_1/eq_{i}_prod_6yb7_pep12.dcd', 2000))
    simulation.step(250000)
    total_sim_time += 250000 * 0.002

    if i == (rounds - 2):
        restraint_atoms = []
        restraint_force = 0
    else:
        restraint_force = restraint_force / 2
        if restraint_force < 1:
            restraint_atoms = alpha_atoms

    positions = simulation.context.getState(getPositions=True).getPositions()

print(total_sim_time)

### Save out PDB #####

integrator = LangevinIntegrator(298 * kelvin, 1 / picosecond, 0.002 * picoseconds)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(positions)
simulation.reporters.append(PDBReporter('binding_1_production_md_500ns_1/output_6yb7_pep12.pdb', 1))
simulation.step(1)

### Production MD #####

restraint_atoms = bbone_atoms
restraint_force = 5
system = make_new_system(prmtop, restraint_atoms, restraint_force, positions)
integrator = LangevinIntegrator(298 * kelvin, 1 / picosecond, 0.002 * picoseconds)
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(positions)
system_xml = XmlSerializer.serialize(system)
with open('binding_1_production_md_500ns_1/output_6yb7_pep12.xml', 'w') as f:
    f.write(system_xml)
simulation.reporters.append(app.StateDataReporter(stdout, 10000, step=True,
                                                  potentialEnergy=True, temperature=True, progress=True,
                                                  remainingTime=True,
                                                  speed=True, totalSteps=1000, separator='\t'))
simulation.reporters.append(DCDReporter('binding_1_production_md_500ns_1/prod_6yb7_pep12.dcd', 10000))
simulation.step(250000000)
total_sim_time += 250000000 * 0.002

print(total_sim_time, 'ns')

