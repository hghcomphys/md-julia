module MDSimulation

include("MDUnits.jl")
include("MDUtils.jl")
include("MDAtoms.jl")
include("MDPotential.jl")
include("MDIntegrator.jl")
include("MDThermostat.jl")
include("MDSystem.jl")
include("MDParams.jl")


using Reexport
@reexport using
	.MDUnits,
	.MDUtils,
	.MDAtoms,
	.MDPotential,
	.MDIntegrator,
	.MDThermostat,
	.MDSystem,
	.MDParams

export simulate!

function simulate_one_step!(system)
	update!(system.integrator, system)
	apply_coupling!(system.thermostat, system)
end

function simulate!(
	system,
	num_steps::Integer = 1,
	output_freq::Integer = 100,
)
	system.atoms.forces = calculate_forces(system.potential, system.atoms)
	print_physical_params(system)
	for step in 1:num_steps
		simulate_one_step!(system)
		if rem(step, output_freq) == 0
			print_physical_params(system)
		end
	end
end

end
