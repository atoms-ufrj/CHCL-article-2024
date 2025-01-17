Github repository associated with the manuscript "Assessing Melting and Solid-Solid Transition Properties of Choline Chloride via Molecular Dynamics Simulations".

The repository is organized into three main folders:

--> LAMMPS_Modification: contains a customized version of the LAMMPS dihedral_style charmm command, which was applied to couple and decouple Lennard-Jones interactions within the CHARMM force field.

--> Unit_Cell: contains the atomic positions of the unit cell of ChCl form beta.

--> Melting_OPLS_DES_0.8: contains the PSCP free energy calculations considering, as an example, the OPLS-DES force field with a charge scaling factor of 0.8. The subfolders are described as: npt - relative free energy as a function of temperature; step0 - parameters of the Gaussian potential; step1 - S->DWF, step2 - DWF->WF, step3 - WF->L.
