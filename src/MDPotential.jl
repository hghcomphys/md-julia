module MDPotential

using ..MDAtoms
using ..MDUtils

export Potential, LJPotential
export calculate_energy, calculate_forces

abstract type Potential{T <: AbstractFloat} end

struct LJPotential{T} <: Potential{T}
	σ::T
	ϵ::T
	r_cutoff::T
end

function calculate_pair_energies(potential::LJPotential, distances)
	λ = potential.σ ./ distances
	λ⁶ = λ .^ 6
	@. 4.0 * potential.ϵ * λ⁶ * (λ⁶ - 1.0)
end

function calculate_energy(potential::LJPotential, atoms)
	natoms = get_natoms(atoms)
	energy = 0.0
	for atom_index in 1:natoms
		distances, _ = calculate_distances(atoms, atom_index)
		energy_per_atom = sum(
			ifelse.(
				calculate_cutoff_masks(distances, potential.r_cutoff),
				calculate_pair_energies(potential, distances),
				0.0,
			),
		)
		# multi-thread with atomic add 
		energy += energy_per_atom
	end
	0.5 * energy
end

function calculate_pair_forces(potential::LJPotential, distances, position_differences)
	λ = potential.σ ./ distances
	λ⁶ = λ .^ 6
	force_factor = @. -24.0 * potential.ϵ / (distances * distances) * λ⁶ * (2.0 * λ⁶ - 1.0)
	add_dim(force_factor) .* position_differences
end

function calculate_forces(potential::LJPotential, atoms)
	natoms = get_natoms(atoms)
	forces = zeros_like(atoms.positions)
	Threads.@threads for atom_index in 1:natoms
		distances, position_differences = calculate_distances(atoms, atom_index)
		pair_forces = ifelse.(
			calculate_cutoff_masks(distances, potential.r_cutoff),
			calculate_pair_forces(potential, distances, position_differences),
			zeros_like(position_differences),
		)
		forces[atom_index, :] = sum(pair_forces, dims = 1)
	end
	return forces
end

end
