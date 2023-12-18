module MDSimulation

include("MDUnits.jl")
include("MDUtils.jl")
include("MDAtoms.jl")
include("MDPotential.jl")
include("MDThermostat.jl")
include("MDIntegrator.jl")
include("MDSystem.jl")
include("MDParams.jl")
include("MDIO.jl")

using Reexport
@reexport using
	.MDUnits,
	.MDUtils,
	.MDAtoms,
	.MDPotential,
	.MDIntegrator,
	.MDThermostat,
	.MDSystem,
	.MDParams,
	.MDIO

export simulate!

function simulate!(
	system,
	num_steps::Integer = 1,
	output_freq::Integer = 100,
	filename = nothing,
)
	system.atoms.forces = calculate_forces(system.potential, system.atoms)
	if !isnothing(filename)
		fio = open(filename, "w")
		dump_xyz(fio, system)
	end
	print_physical_params(system)
	for step in 1:num_steps
		simulate_one_step!(system.integrator, system)
		if step % output_freq == 0
			print_physical_params(system)
			if !isnothing(filename)
				dump_xyz(fio, system)
			end
		end
	end
	if !isnothing(filename)
		close(fio)
	end
end

end
