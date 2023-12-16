module MDIntegrator

using ..MDAtoms
using ..MDPotential

export
	Integrator,
	VelocityVerlet,
	update!

abstract type Integrator{T <: AbstractFloat} end

mutable struct VelocityVerlet{T} <: Integrator{T}
	time_step::T
	step::Integer
	time::T
end

VelocityVerlet{T}(
	time_step,
) where {T <: AbstractFloat} = VelocityVerlet(T(time_step), 0, T(0.0))

function verlet_new_positions(positions, velocities, forces, time_step)
	@. positions + (velocities + 0.5 * forces * time_step) * time_step
end

function verlet_new_velocities(velocities, forces, new_forces, time_step)
	@. velocities + 0.5 * (forces + new_forces) * time_step
end

function update!(::VelocityVerlet, system)
	atoms = system.atoms
	time_step = system.integrator.time_step

	new_positions = verlet_new_positions(
		atoms.positions, atoms.velocities, atoms.forces, time_step,
	)
	# system.atoms.positions = new_positions
	system.atoms.positions = shift_inside_box(new_positions, atoms.lattice)

	new_forces = calculate_forces(system.potential, system.atoms)
	new_velocities = verlet_new_velocities(
		atoms.velocities, atoms.forces, new_forces, time_step,
	)
	system.atoms.forces = new_forces
	system.atoms.velocities = new_velocities

	system.integrator.step += 1
	system.integrator.time += time_step
end

end
