module MDSystem

using ..MDAtoms
using ..MDPotential
using ..MDIntegrator
using ..MDThermostat

export System

mutable struct System{T <: AbstractFloat}
	atoms::Atoms{T}
	potential::Potential{T}
	integrator::Integrator{T}
end

end
