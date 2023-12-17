module MDParams

using Printf
using ..MDAtoms
using ..MDUnits
using ..MDUtils
using ..MDPotential

export
	calculate_potential_energy,
	calculate_kinetic_energy,
	calculate_total_energy,
	calculate_temperature,
	calculate_pressure,
	print_physical_params

function calculate_potential_energy(system)
	return calculate_energy(system.potential, system.atoms)
end

function calculate_kinetic_energy(system)
	atoms = system.atoms
	return sum(@. 0.5 * (atoms.masses * atoms.velocities * atoms.velocities))
end

function calculate_total_energy(system)
	return (
		calculate_potential_energy(system)
		+
		calculate_kinetic_energy(system)
	)
end

function calculate_temperature(system)
	temperature_kernel(system.atoms.velocities, system.atoms.masses)
end

function calculate_pressure(system)
	atoms = system.atoms
	virial = (
		2 * calculate_kinetic_energy(system)
		+
		sum(atoms.positions .* atoms.forces)
	)
	return virial / (3 * get_volume(atoms))
end

function calculate_center_mass_velocity(system)
	calculate_center_of_mass(system.atoms.velocities, system.atoms.masses)
end

function calculate_center_mass_position(system)
	calculate_center_of_mass(system.atoms.positions, system.atoms.masses)
end

function print_physical_params(system)
	@printf(
		"%-10d t[ps]:%-10.5f Temp[K]:%-15.10f Etot[Ha]:%-15.10f Epot[Ha]:%-15.10f Pres[kb]:%-15.10f\n",
		system.integrator.step,
		system.integrator.time * TO_PS,
		calculate_temperature(system),
		calculate_total_energy(system),
		calculate_potential_energy(system),
		calculate_pressure(system) * TO_KB,
	)
	# println("Vcom: ", calculate_center_mass_velocity(system))
	# println("Rcom: ", calculate_center_mass_position(system))
	# println("Drifting force: ", sum(system.atoms.forces, dims = 1))
end

end
