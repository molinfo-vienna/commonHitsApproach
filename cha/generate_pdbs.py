"""
Created on Mon Oct 25 08:00:03 2016
@author: Marcus Wieder as mwieder
department of pharmaceutical chemistry
university vienna

"""

# =============================================================================================
# IMPORTS
# =============================================================================================

import MDAnalysis

def generate_pdbs_from_charmm_input(trajfile, psffile, pdb_output_directory, ligand_3_letter_code):
    print( "Generating pdb files from:")
    print( trajfile)
    print( psffile)
    print( ligand_3_letter_code)
    print( "in " + str(pdb_output_directory))
    u = MDAnalysis.Universe(psffile,  trajfile)
    protein_ligand = u.select_atoms("resname " + str(ligand_3_letter_code) + " or protein")

    print( 'Reading trajectory and generating pdbs ... ')
    count = 0
    for ts in u.trajectory:
        print( ts)
        count += 1
        protein_ligand.write(str(pdb_output_directory) + '/lp_' + str(count) + '.pdb')

