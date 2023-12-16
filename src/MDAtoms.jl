module MDAtoms

using LinearAlgebra
using FieldProperties
using ..MDUtils

export
	Atoms,
	calculate_distances,
	calculate_cutoff_masks,
	get_volume,
	get_natoms,
	shift_inside_box

mutable struct Atoms{T <: AbstractFloat}
	lattice::Array{T, 2} # SArray?
	positions::Array{T, 2}
	velocities::Array{T, 2}
	forces::Array{T, 2}
	masses::Array{T, 1}
end

function get_natoms(atoms::Atoms)
	size(atoms.positions, 1)
end

function get_volume(atoms::Atoms)
	prod(diag(atoms.lattice))
end

function apply_pbc!(position_differences, lattice)
	for i in 1:3
		dx = position_differences[:, i]
		l = lattice[i, i]
		dx = ifelse.(dx .> 0.5 * l, dx .- l, dx)
		dx = ifelse.(dx .< -0.5 * l, dx .+ l, dx)
		@assert sum(abs.(dx) .> l) == 0
		position_differences[:, i] = dx
	end
end

function calculate_cutoff_masks(distances, r_cutoff)
	(distances .<= r_cutoff) .& (distances .!= 0.0)
end

function calculate_distances(atoms::Atoms, atom_index::Integer)
	atom_position = reshape(atoms.positions[atom_index, :], (1, 3))
	position_differences = atom_position .- atoms.positions
	apply_pbc!(position_differences, atoms.lattice)
	distances = sqrt.(sum(position_differences .^ 2, dims = 2))[:, 1]
	return distances, position_differences
end

function shift_inside_box(positions, lattice)
	positions .% reshape(diag(lattice), (1, 3))
end

end
