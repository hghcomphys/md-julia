module MDThermostat

using ..MDAtoms
using ..MDUtils

export
	Thermostat,
	BrendsenThermostat,
	apply_coupling!

abstract type Thermostat{T <: AbstractFloat} end

struct BrendsenThermostat{T} <: Thermostat{T}
	target_temperature::T
	time_constant::T
end

function apply_coupling!(thermostat::BrendsenThermostat, system)
	scaling_factor = 1.0 / sqrt(
		1.0
		+
		(system.integrator.time_step / thermostat.time_constant)
		*
		(temperature_kernel(system.atoms.velocities, system.atoms.masses) / thermostat.target_temperature - 1.0),
	)
	system.atoms.velocities .*= scaling_factor
end

end
