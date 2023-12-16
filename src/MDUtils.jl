module MDUtils

using Random
using ..MDUnits

export
	add_dim,
	zeros_like,
	temperature_kernel,
	generate_random_velocities,
	calculate_center_of_mass

add_dim(x::Array) = reshape(x, (size(x)..., 1))

zeros_like(x::Array) = zeros(eltype(x), size(x))

function temperature_kernel(velocities, masses)
	natoms = size(masses, 1)
	return (
		sum(@. (masses * velocities * velocities))
		/
		(3 * natoms * BOLTZMANN_CONSTANT)
	)
end

function calculate_center_of_mass(array, masses)
	sum(add_dim(masses) .* array, dims = 1) / sum(masses)
end

function generate_random_velocities(temperature, masses, seed = 2023)
	"""Generate Maxwell-Boltzmann distributed random velocities."""
	Random.seed!(seed)
	natoms = size(masses, 1)
	velocities = randn((natoms, 3))
	velocities *= sqrt.(
		temperature ./ temperature_kernel(velocities, masses),
	)
	velocities .-= calculate_center_of_mass(velocities, masses)
	return velocities
end

end
