module MDIO

export dump_xyz

using Printf
using ..MDAtoms
using ..MDUnits

function dump_xyz(fio::IOStream, system)
	atoms = system.atoms
	natoms = get_natoms(atoms)
	@printf(fio, "%d\n\n", natoms)	
	for atom_index in 1:natoms
		@printf(
			fio, 
			"X %1.10f %1.10f %1.10f\n", 
			atoms.positions[atom_index, 1] * TO_ANG,
			atoms.positions[atom_index, 2] * TO_ANG,
			atoms.positions[atom_index, 3] * TO_ANG
		)
	end
end

end
