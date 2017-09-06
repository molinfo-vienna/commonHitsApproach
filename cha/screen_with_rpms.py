import re, glob, os
import argparse
import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Pharm as Pharm
import CDPL.Biomol as Biomol
from collections import defaultdict
import subprocess
import time


ftype_names = { Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA', Pharm.FeatureType.H_BOND_DONOR : 'HBD', Pharm.FeatureType.POS_IONIZABLE : 'PI', Pharm.FeatureType.NEG_IONIZABLE : 'NI', Pharm.FeatureType.AROMATIC : 'AR', Pharm.FeatureType.HYDROPHOBIC : 'H', Pharm.FeatureType.X_VOLUME : 'XV'  }


class Screen():
    def __init__(self, rpm_directory, database):
        self.rpm_directory = rpm_directory
        self.database = database


    def generate_screening_database(self):

            def setupMolecule(mol):
                Chem.perceiveComponents(mol, False)
                Chem.perceiveSSSR(mol, False)
                Chem.setRingFlags(mol, False)
                Chem.calcImplicitHydrogenCounts(mol, False)
                Chem.perceiveHybridizationStates(mol, False)
                Chem.setAromaticityFlags(mol, False)
                Chem.calcCIPPriorities(mol, False)
                Chem.calcAtomCIPConfigurations(mol, False)
                Chem.calcBondCIPConfigurations(mol, False)

            def process(sdf_file, psd_file_path):

                ifs = Base.FileIOStream(sdf_file, 'r')

                reader = Chem.SDFMoleculeReader(ifs)

                mol = Chem.BasicMolecule()

                Chem.setMultiConfImportParameter(reader, True)

                psd_creator = Pharm.PSDScreeningDBCreator(psd_file_path, Pharm.PSDScreeningDBCreator.CREATE, True)
                i = 0;
                t0 = time.clock()

                while reader.read(mol):
                    setupMolecule(mol)

                    psd_creator.process(mol)
                    i += 1

                    if i % 100 == 0:
                        print 'Processed ' + str(i) + ' molecules (' + str(time.clock() - t0), 's elapsed)...'
                        t0 = time.clock()

                    mol.clear()

                print ''
                print '-- Summary --'
                print 'Molecules processed: ' + str(psd_creator.numProcessed)
                print 'Molecules rejected: ' + str(psd_creator.numRejected)
                print 'Molecules deleted: ' + str(psd_creator.numDeleted)
                print 'Molecules inserted: ' + str(psd_creator.numInserted)

                psd_creator.close()



    def read_rpms(args):
        def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
            return [int(text) if text.isdigit() else text.lower()
                for text in re.split(_nsre, s)]

        def read_in_ph(ph_path, output_dir_path):

            fr = Pharm.PMLPharmacophoreReader(Base.FileIOStream(ph_path))
            ph = Pharm.BasicPharmacophore()
            fr.read(ph)
            ph.pml_path = ph_path
            ph.dir_path = output_dir_path

            return ph


        dir_list = glob.glob(str(args.input) + '*/')
        rpm_map = defaultdict(list)

        for dir_path in sorted(dir_list, key=natural_sort_key):


            a,b = os.path.split(str(dir_path))
            a,fv = os.path.split(str(a))
            if str(fv) == 'ligand':
                continue
            #print dir_path

            pml_list = glob.glob(str(dir_path) + '/*pml')

            for pml in sorted(pml_list, key=natural_sort_key):
                a,b = os.path.split(str(pml))
                basename = b.replace('.pml', '').replace('ph_','')
                ph = read_in_ph(pml, dir_path)
                ph.basename = basename
                rpm_map[fv].append(ph)

        return rpm_map

    def filter_rpms(rpm_map):

        # here we need to decide on the filtering criteria
        screening_ph = rpm_map
        return screening_ph

    def start_screening(screening_ph, args):
        print 'Removing old files ...'
        for ph_list in screening_ph.values():
            for ph in ph_list:
                print 'Submitting ...'
                try:
                    os.remove(str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.log')
                    os.remove(str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.sdf')
                except OSError:
                    pass

                print 'qsub -o ' + str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.log  -N psd_screen_' + str(ph.basename), '/data/shared/projects/scripts/sge_cdpl_screening.sh', '-q', str(ph.pml_path), '-d', args.database , '-o', str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.sdf'


                output = subprocess.Popen(['qsub', '-o', str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.log',  '-N', 'psd_screen_' + str(ph.basename), '/data/shared/projects/scripts/sge_cdpl_screening.sh', '-q', str(ph.pml_path), '-d', args.database , '-o', str(ph.dir_path) + '/hitlist_' + str(ph.basename) + '.sdf'], stdout=subprocess.PIPE, shell=False)

            for line in output.stdout:
                print line.rstrip()


