import sys
import argparse
import MDAnalysis
import os
import pandas as pd
import numpy as np
import random
sys.path.insert(0,'.' )
import cha

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest="trajectory", required=True, help="input trajectory (dcd) file")
    parser.add_argument('-p', dest="psf", required=True, help="input psf file")
    parser.add_argument('-l', dest="ligand_three_letter_code" ,required=True, help="ligand 3 letter code")
    parser.add_argument('-o', dest='output', required=True, help='specify filname for representative pharmacophore models. Will be saved as csv. Do not add ".csv" to the name.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_arguments()
    trajfile = args.trajectory

    df_constructor = dict()
    psffile = args.psf
    u = MDAnalysis.Universe(psffile, trajfile)
    protein_ligand = u.select_atoms("resname " + str(args.ligand_three_letter_code) + " or protein")
    length = 1
    interaction_at_all_ts = dict()
    pdb_list = list()
    for ts in u.trajectory:
        rand = random.random()
        print( '***************************')
        print( str(length))
        name = './proteins/protein_ligand_'+str(rand)+'.pdb'
        protein_ligand.write(name)
        pdb_list.append(name)
        length += 1
        if length == 6:
            break
    map_of_pdbs = {}
    map_of_pdbs['frame'] = pdb_list

    rpms = cha.generate_rpms.generate_rpms_factory(map_of_pdbs, args.ligand_three_letter_code)
    ph_map,unique_str = rpms.read_pdb_and_generate_ph()
    rpm_map = rpms.generate_rpms(ph_map)
    df = pd.DataFrame(rpm_map.items())
    print("df",df)
    df.to_csv(args.output+".csv")
    text_file = open(args.output+"_unique_pha_vector.txt", "w")
    text_file.write(unique_str)
    text_file.close()

    for file in pdb_list:
        os.remove(file)

