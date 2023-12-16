include("../src/MDSimulation.jl")

using ExtXYZ: read_frame
using .MDSimulation

const T = Float64

println("Molecular Dynamics Simulation in Julia")
println("Example: Lennard-Jones Potential")
println("Precision: ", T)

mass = 4.003 * FROM_AMU #  
temperature = 300.0  # Kelvin
σ = 2.5238 * FROM_ANG  # Bohr
ϵ = 4.7093e-04 * FROM_EV  # Hartree
r_cutoff = 6.3095 * FROM_ANG # Angstrom
time_step = 0.5 * FROM_FS  # hbar/Hartree

data = read_frame("He100.xyz")
lattice = FROM_ANG * data["cell"]
positions = shift_inside_box(FROM_ANG * transpose(data["arrays"]["pos"]), lattice)
natoms = data["N_atoms"]
masses = mass * ones(T, (natoms,))
velocities = generate_random_velocities(temperature, masses)
forces = similar(velocities)
atoms = Atoms{T}(lattice, positions, velocities, forces, masses)
println("Number of atoms: ", get_natoms(atoms))
# println("distances: ", calculate_distances(atoms, 2))

potential = LJPotential{T}(σ, ϵ, r_cutoff)
# println("potential: ", potential)
# println("energy: ", calculate_energy(potential, atoms))
# println("forces: ", calculate_forces(potential, atoms))

integrator = VelocityVerlet{T}(time_step)
thermostat = BrendsenThermostat{T}(temperature, 100 * time_step)
system = System{T}(atoms, potential, integrator, thermostat)

simulate!(system, 10_000, 100, "dump.xyz")

