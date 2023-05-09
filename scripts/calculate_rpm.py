import sys, re
import argparse
import math
import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import CDPL.Biomol as Biomol
import CDPL.MolProp as MolProp

import MDAnalysis
import pandas as pd
import numpy as np
import matplotlib.pylab as pl
#import seaborn as sbn
from collections import defaultdict
import random
sys.path.insert(0,'.' )
import cha

ftype_names = { Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA', Pharm.FeatureType.H_BOND_DONOR : 'HBD', Pharm.FeatureType.POSITIVE_IONIZABLE : 'PI', Pharm.FeatureType.NEGATIVE_IONIZABLE : 'NI', Pharm.FeatureType.AROMATIC : 'AR', Pharm.FeatureType.HYDROPHOBIC : 'H'}#, Pharm.FeatureType.X_VOLUME : 'XV'  }


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
        for text in re.split(_nsre, s)]

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', dest="trajectory", required=True, help="input trajectory (dcd) file")
    parser.add_argument('-p', dest="psf", required=True, help="input psf file")
    parser.add_argument('-l', dest="ligand_three_letter_code" ,required=True, help="ligand 3 letter code")
    parser.add_argument('-o', dest='output', required=True, help='specify filname for representative pharmacophore models')
    args = parser.parse_args()
    return args


def generate_key(ftr):
    key_list = []
    first_atom = Pharm.getSubstructure(ftr).atoms[0]
    base = str(ftype_names[Pharm.getType(ftr)]) + '[' + str(Biomol.getResidueCode(first_atom)) + '_' + str(Biomol.getResidueSequenceNumber(first_atom)) + '_' + str(Biomol.getChainID(first_atom))
    atoms_list = []
    for a in Pharm.getSubstructure(ftr).atoms:
        if Biomol.hasSerialNumber(a) == False:
            continue

        atom_id = str(Biomol.getSerialNumber(a))
        atom_id += ':'
        atom_id += str(Chem.getSymbol(a))
        atoms_list.append(atom_id)

    atom_key = ""
    for k in sorted(atoms_list, key=natural_sort_key, reverse=True):
        atom_key += '_' + k

    key = base + atom_key + ']'
    #print key
    return key

def outputInteractions(lig_pharm, env_pharm, interactions, df_constructor):
    i = 0

    interaction_at_ts = dict()

    for lig_ftr in lig_pharm:
        if Pharm.hasSubstructure(lig_ftr) == False:
           continue
        elif ftype_names[Pharm.getType(lig_ftr)] == 'XV':
           continue
        elif len(interactions.getValues(lig_ftr)) < 1:
           continue
        ligand_key = generate_key(lig_ftr)
        print( 'Ligand feature : ' + str(ligand_key) + ' interacts with: ')

        env_ftrs = interactions.getValues(lig_ftr)
        if ligand_key in df_constructor:
            dic_of_env_key = df_constructor[ligand_key]
        else:
            dic_of_env_key = {}

        dic_of_env_key_at_ts = {}
        for env_ftr in env_ftrs:
            if Pharm.hasSubstructure(env_ftr) == False:
               continue
            elif ftype_names[Pharm.getType(lig_ftr)] == 'XV':
               continue
            env_key = generate_key(env_ftr)
            if env_key in dic_of_env_key:
                dic_of_env_key_at_ts[env_key] = 1
                dic_of_env_key[env_key] += 1
            else:
                dic_of_env_key[env_key] = 1
                dic_of_env_key_at_ts[env_key] = 1


            print( ' - ' + str(env_key))

        df_constructor[ligand_key] = dic_of_env_key
        interaction_at_ts[ligand_key] = dic_of_env_key_at_ts


    return df_constructor, interaction_at_ts


def generate_ph(pdb, args, df_constructor, ts):

    ifs = Base.FileIOStream(pdb, 'r')
    tlc = args.ligand_three_letter_code
    pdb_reader = Biomol.PDBMoleculeReader(ifs)
    pdb_mol = Chem.BasicMolecule()

    print( '- Reading input: ', pdb, ' ...')

    if not pdb_reader.read(pdb_mol):
        print( '!! Could not read input molecule')
        return

    print( '- Processing macromolecule', pdb, ' ...')

    i = 0

    while i < pdb_mol.getNumBonds():
        bond = pdb_mol.getBond(i)

        if MolProp.isMetal(bond.atoms[0]) or MolProp.isMetal(bond.atoms[1]):
            pdb_mol.removeBond(i)
        else:
            i += 1

    for a in pdb_mol.atoms:
        Chem.setImplicitHydrogenCount(a, 0)

    Chem.calcImplicitHydrogenCounts(pdb_mol, True)
    Chem.perceiveHybridizationStates(pdb_mol, True)
    Chem.makeHydrogenComplete(pdb_mol)
    Chem.setAtomSymbolsFromTypes(pdb_mol, False)
    Chem.calcImplicitHydrogenCounts(pdb_mol, True)
    Biomol.setHydrogenResidueSequenceInfo(pdb_mol, False)
    Chem.setRingFlags(pdb_mol, True)
    Chem.setAromaticityFlags(pdb_mol, True)
    Chem.calcHydrogen3DCoordinates(pdb_mol, True)
    Chem.calcFormalCharges(pdb_mol, True)
    Pharm.prepareForPharmacophoreGeneration(pdb_mol)  
    ligand = Chem.Fragment()

    print( '- Extracting ligand ', tlc, ' ...')

    for atom in pdb_mol.atoms:
        if Biomol.getResidueCode(atom) == tlc:
            Biomol.extractResidueSubstructure(atom, pdb_mol, ligand, False)
            break

    if ligand.numAtoms == 0:
        print( '!! Could not find ligand', tlc, 'in input file')
        return

    Chem.perceiveSSSR(ligand, True)

    lig_env = Chem.Fragment()

    Biomol.extractEnvironmentResidues(ligand, pdb_mol, lig_env, 7.0)
    Chem.perceiveSSSR(lig_env, True)
    print( '- Constructing pharmacophore ...')
    
    #### new
    lig_pharm = Pharm.BasicPharmacophore()
    env_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.DefaultPharmacophoreGenerator()

    # Pharm.prepareForPharmacophoreGeneration(ligand)
    pharm_gen.generate(ligand, lig_pharm)
    # Pharm.prepareForPharmacophoreGeneration(lig_env)
    pharm_gen.generate(lig_env, env_pharm)
    #Pharm.FilePMLFeatureContainerWriter('./test/lig_ph_' + str(ts) + '.pml').write(lig_pharm)

    analyzer = Pharm.DefaultInteractionAnalyzer()
    interactions = Pharm.FeatureMapping()
    analyzer.analyze(lig_pharm, env_pharm, interactions)
    df_constructor, interaction_at_ts = outputInteractions(lig_pharm, env_pharm, interactions, df_constructor)
	#Chem.FileSDFMolecularGraphWriter('./test/ligand_' + str(ts) + '.sdf').write(ligand)


    return df_constructor, interaction_at_ts
    
    # #### old
    # lig_pharm = Pharm.BasicPharmacophore()
    # env_pharm = Pharm.BasicPharmacophore()
    # pharm_gen = Pharm.DefaultPharmacophoreGenerator(True)
    # pharm_gen.generate(ligand, lig_pharm)
    # pharm_gen.generate(lig_env, env_pharm)
    # #Pharm.FilePMLFeatureContainerWriter('./test/lig_ph_' + str(ts) + '.pml').write(lig_pharm)

    # analyzer = Pharm.DefaultInteractionAnalyzer()
    # interactions = Pharm.FeatureMapping()
    # analyzer.analyze(lig_pharm, env_pharm, interactions)
    # df_constructor, interaction_at_ts = outputInteractions(lig_pharm, env_pharm, interactions, df_constructor)
	# #Chem.FileSDFMolecularGraphWriter('./test/ligand_' + str(ts) + '.sdf').write(ligand)


    # return df_constructor, interaction_at_ts

def plot_frequency_of_features(args, interaction_at_all_ts, set_of_ligand_interaction_partner, set_of_protein_interaction_partner, length, binning):

    full_interactions = defaultdict(dict)

    for ligand_interaction in sorted(set_of_ligand_interaction_partner):

        for protein_interaction in sorted(set_of_protein_interaction_partner):
            a = np.zeros(length)
            full_interactions[ligand_interaction][protein_interaction] = a

    for ts in interaction_at_all_ts.keys():

        for ligand_interaction in sorted(set_of_ligand_interaction_partner):

            for protein_interaction in sorted(set_of_protein_interaction_partner):

                if ligand_interaction in interaction_at_all_ts[ts]:
                    if protein_interaction in interaction_at_all_ts[ts][ligand_interaction]:
                        full_interactions[ligand_interaction][protein_interaction][ts]=1

    number_of_plots = 1
    # print("full_interactions",full_interactions)
    for ligand_interaction in sorted(set_of_ligand_interaction_partner):

        for protein_interaction in sorted(set_of_protein_interaction_partner):
            # print( full_interactions[ligand_interaction][protein_interaction])
            check = np.sum(full_interactions[ligand_interaction][protein_interaction])
            if check != 0:
                number_of_plots += 1

    axes_list = []
    f, (axes_list) = pl.subplots(number_of_plots, sharex=True, sharey=False, figsize=(18, 2*number_of_plots))
    axes_list = list(axes_list.tolist())
    for ligand_interaction in sorted(set_of_ligand_interaction_partner):
        for protein_interaction in sorted(set_of_protein_interaction_partner):
            # print( full_interactions[ligand_interaction][protein_interaction])
            check = np.sum(full_interactions[ligand_interaction][protein_interaction])
            if check != 0:
                # print( str(ligand_interaction) + ' - ' + str(protein_interaction) + ': ' + str(check))
                binned_ts_list = []
                for j in range(0, len(full_interactions[ligand_interaction][protein_interaction]),binning):
                    binned_ts_list.append(full_interactions[ligand_interaction][protein_interaction][j:j+binning:1].sum())
                # print( binned_ts_list)
                specific_axis = axes_list.pop()
                x= np.linspace(1, len(binned_ts_list),len(binned_ts_list))
                specific_axis.plot(x,binned_ts_list,  linewidth=2, color='blue')
                specific_axis.set_title('Interaction between :' + ligand_interaction + ' and ' + protein_interaction, fontsize=12)
                specific_axis.set_ylim(0, float(binning) + 0.5)
                specific_axis.set_xlim(1,len(binned_ts_list))
                specific_axis.fill_between(x,binned_ts_list,facecolor='blue', alpha=0.2)

    f.subplots_adjust(hspace=2.5)
    f.tight_layout()
    f.savefig(args.output_2, format='svg')


def generate_interaction_map(set_of_ligand_interaction_partner, set_of_protein_interaction_partner, args):


    container = defaultdict(list)
    # print( interaction_at_all_ts)

    for ligand_interaction in sorted(set_of_ligand_interaction_partner):
        ano = []
        for protein_interaction in sorted(set_of_protein_interaction_partner):
            if protein_interaction in df_constructor[ligand_interaction]:
                ano.append(df_constructor[ligand_interaction][protein_interaction])
            else:
                ano.append(0)
        container[ligand_interaction] = ano
        print("ano",ano)
        print("container",container)

    df = pd.DataFrame(container, index=sorted(list(set_of_protein_interaction_partner)))
    n_rows, m_columns = df.shape
    y_axis=[]
    for i in range(0, n_rows):
        y_axis.append(i)

    x_axis=[]
    for i in range(0, m_columns):
        x_axis.append(i)

    xminor_ticks = np.arange(0.5, len(x_axis),1)
    yminor_ticks = np.arange(0.5, len(y_axis),1)
    figsize = (20,20)

    fig = pl.figure( figsize = figsize )
    x_label = 'Ligand Atoms'
    y_label = 'Protein Atoms'
    titel = 'test'
    ax = fig.add_subplot(121)
    ax.set_ylabel(x_label)
    ax.set_xlabel(y_label)
    ax.set_title(titel)
    ax.get_frame_on()
    ax.set_xticks(x_axis)
    ax.set_xticklabels(df.columns, rotation='vertical' )
    ax.set_yticks(y_axis)
    ax.set_yticklabels(df.index)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(yminor_ticks, minor=True)
    c = ax.imshow(df, interpolation="nearest",origin="lower")
    fig.colorbar(c)
    # pl.grid(b=True, which='minor', color='0.65',linestyle='-')
    pl.savefig(args.output_1, format='svg')



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
        name = './protein_ligand_'+str(rand)+'.pdb'
        protein_ligand.write(name)
        pdb_list.append('./'+name)
        length += 1
        if length == 6:
            break
    map_of_pdbs = {}
    map_of_pdbs['mr1'] = pdb_list

    rpms = cha.generate_rpms.generate_rpms_factory(map_of_pdbs, args.ligand_three_letter_code)
    ph_map = rpms.read_pdb_and_generate_ph()
    rpms.generate_rpms(ph_map)