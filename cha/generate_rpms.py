import  re, os
import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import CDPL.Biomol as Biomol
import CDPL.MolProp as MolProp
from collections import defaultdict


# ftype_names = { Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA', Pharm.FeatureType.H_BOND_DONOR : 'HBD', 
#                Pharm.FeatureType.POS_IONIZABLE : 'PI', Pharm.FeatureType.NEG_IONIZABLE : 'NI', 
#                Pharm.FeatureType.AROMATIC : 'AR', Pharm.FeatureType.HYDROPHOBIC : 'H', Pharm.FeatureType.X_VOLUME : 'XV'  }

ftype_names = { Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA', Pharm.FeatureType.H_BOND_DONOR : 'HBD', 
               Pharm.FeatureType.POSITIVE_IONIZABLE : 'PI', Pharm.FeatureType.NEGATIVE_IONIZABLE : 'NI', 
               Pharm.FeatureType.AROMATIC : 'AR', Pharm.FeatureType.HYDROPHOBIC : 'H', Pharm.FeatureType.EXCLUSION_VOLUME : 'XV'  }



class generate_rpms_factory():
    def __init__(self, map_of_pdbs,  ligand_3_letter_code):
        self.unique_feature_vector = set()
        self.list_of_pdbs = map_of_pdbs # map[mr] = list(pdb1_file_path, pdb2_file_path)
        self.ligand_3_letter_code = ligand_3_letter_code





    def read_pdb_and_generate_ph(self):
        def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
            return [int(text) if text.isdigit() else text.lower()
                for text in re.split(_nsre, s)]


        def generate_key(ftr):
            first_atom = Pharm.getSubstructure(ftr).atoms[0]
            base = str(ftype_names[Pharm.getType(ftr)]) + '[' + str(Biomol.getResidueCode(first_atom)) + '_' + str(Biomol.getResidueSequenceNumber(first_atom)) + '_' + str(Biomol.getChainID(first_atom))
            atoms_list = []
            for a in Pharm.getSubstructure(ftr).atoms:
                if Biomol.hasSerialNumber(a) == False:
                    continue

                atom_id = str(Biomol.getSerialNumber(a))
                atoms_list.append(atom_id)

            atom_key = ""
            for k in sorted(atoms_list, key=natural_sort_key, reverse=True):
                atom_key += '_' + k

            key = base + atom_key + ']'
            return key





        def generate_ph(pdb, key):

            ifs = Base.FileIOStream(pdb, 'r')
            tlc = self.ligand_3_letter_code
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
            lig_pharm = Pharm.BasicPharmacophore()
            env_pharm = Pharm.BasicPharmacophore()
            ### new
            #### new
            pharm_gen = Pharm.DefaultPharmacophoreGenerator()

            # Pharm.prepareForPharmacophoreGeneration(ligand)
            pharm_gen.generate(ligand, lig_pharm)
            # Pharm.prepareForPharmacophoreGeneration(lig_env)
            pharm_gen.generate(lig_env, env_pharm)
            ###
            analyzer = Pharm.DefaultInteractionAnalyzer()
            interactions = Pharm.FeatureMapping()
            analyzer.analyze(lig_pharm, env_pharm, interactions)

        #------------------------- XVOLS

            int_env_ftrs = Pharm.FeatureSet()
            Pharm.getFeatures(int_env_ftrs, interactions, False)
            int_core_ftrs = Pharm.FeatureSet()
            Pharm.getFeatures(int_core_ftrs, interactions, True)
            int_pharm = Pharm.BasicPharmacophore(int_core_ftrs)

            for ftr in int_env_ftrs:
                if Pharm.getType(ftr) == Pharm.FeatureType.H_BOND_DONOR or Pharm.getType(ftr) == Pharm.FeatureType.H_BOND_ACCEPTOR:
                    Pharm.setTolerance(ftr, 1.0)
                else:
                    Pharm.setTolerance(ftr, 1.5)

            Pharm.createExclusionVolumes(int_pharm, int_env_ftrs, 0.0, 0.1, False)
            int_env_ftr_atoms = Chem.Fragment()
            Pharm.getFeatureAtoms(int_env_ftrs, int_env_ftr_atoms)
            int_residue_atoms = Chem.Fragment()
            Biomol.extractResidueSubstructures(int_env_ftr_atoms, lig_env, int_residue_atoms, True)
            Chem.makeHydrogenDeplete(int_residue_atoms)

            def isAlphaAtom(atom):
                return Biomol.getResidueAtomName(atom) == 'CA'

            Chem.removeAtomsIfNot(int_residue_atoms, isAlphaAtom)
            Pharm.createExclusionVolumes(int_pharm, int_residue_atoms, Chem.Atom3DCoordinatesFunctor(), 1.0, 2.0, False)

            features_in_ph = []
            for int_ftr in int_pharm:
                if Pharm.hasSubstructure(int_ftr) == False:
                       continue
                elif ftype_names[Pharm.getType(int_ftr)] == 'XV':
                       continue
                feature_id = generate_key(int_ftr)
                features_in_ph.append(str(feature_id))
                self.unique_feature_vector.add(str(feature_id))

            int_pharm.fv = features_in_ph
            int_pharm.path_to_pdb = pdb

            return int_pharm




        ###################################################

        ph_map = {}

        for mr  in sorted(self.list_of_pdbs.keys(), key=natural_sort_key):
            count = 0
            pdb_list = self.list_of_pdbs[mr]
            for pdb in pdb_list:
                count += 1

                key = str(mr) + '_' + str(count)
                ph_map[key] = generate_ph(pdb,  key)
                print( '**********************************')

        print( 'Unique feature vector: ' + str(self.unique_feature_vector))
        return ph_map,str(self.unique_feature_vector)



    def generate_rpms(self, ph_map):
        def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
            return [int(text) if text.isdigit() else text.lower()
                for text in re.split(_nsre, s)]


        rpm_map = defaultdict(list)
        unique_feature_list = list(sorted(self.unique_feature_vector, key=natural_sort_key))
        for key in sorted(ph_map, key=natural_sort_key):
            feature_vector = []

            for unique_feature in unique_feature_list:
                if unique_feature in ph_map[key].fv:
                   feature_vector.append(1)
                else:
                    feature_vector.append(0)

            fv = ''.join(str(feature_vector)).replace('[','').replace(']','').replace(',','').replace(' ','')
            print( fv)
            rpm_map[fv].append(key)
        return rpm_map


    def write_ph_for_rpms(self, rpm_maps, output_directory):
        for fv in rpm_maps:
            directory = output_directory + '/' + str(fv)
            if not os.path.exists(directory):
                os.makedirs(directory)
            for ph_key in rpm_maps[fv]:
                ph_to_write =  directory + '/ph_' +str(fv) + '_' + str(ph_key) + '.pml'
                print( '- Writing pharmacophore: ' + str(ph_to_write))
                Pharm.PMLFeatureContainerWriter(Base.FileIOStream(ph_to_write, 'w')).write(rpm_maps[ph_key])



