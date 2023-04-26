#16 November 2016

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import os.path
from rdkit import Chem
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate and save the AUC values for a given VS hitlist. Allow to draw the ROC.")
    parser.add_argument('-i', dest="input", required=True, help="Input sdf files (Hitlist)")
    parser.add_argument('-o', dest="output", required=False, help="Output ROC-CURVE as svg file")
    parser.add_argument('-p', dest="activity_property_name", required=True, help="Name of the activity property")
    parser.add_argument('-na', dest="nb_act", required=True, help="Number of actives in the screened database", type=int)
    parser.add_argument('-nd', dest="nb_dec", required=True, help="Number of decoys in the screened database", type=int)
    parser.add_argument('-v', dest="verbose", help="print additional information", action="store_true", default=False)
    args = parser.parse_args()
    return args


def reading_input(args):
    molecules = Chem.SDMolSupplier(args.input, removeHs=False, sanitize=False)
    return molecules

def calculating_TP_and_FP(molecules, args):

    active = 0
    decoy = 0
    true_positive_rate_for_every_mol = [0.0]
    false_positive_rate_for_every_mol = [0.0]
    if args.verbose == True:
        print( 'Len of hitlist: ' + str(len(molecules)))
    if len(molecules) < 1:
	    print( '0, 0, 0, 0, 0, 0, 0, 0')
    else:
        for molecule in sorted(molecules, key=lambda x: x.GetProp('Score'), reverse=True):
            if float(molecule.GetProp('Score').strip()) == float(0.0):
                print( 'Score Zero!')
                print( molecule.GetProp('Score').strip())
            if molecule.GetProp(args.activity_property_name) == 'active':
                if args.verbose == True:
                    print( 'Active at position: ' + str(active+decoy) + ' in hitlist ...')
                active += 1
            elif molecule.GetProp(args.activity_property_name) == 'decoy':
                decoy += 1
            true_positive_rate = float(active) / float(args.nb_act)
            false_positive_rate = float(decoy) / float(args.nb_dec)
            true_positive_rate_for_every_mol.append(true_positive_rate)
            false_positive_rate_for_every_mol.append(false_positive_rate)

        false_positive_rate_for_every_mol.append(1)
        true_positive_rate_for_every_mol.append(1)
        precision = true_positive_rate/ (true_positive_rate + false_positive_rate)
        accuracy = float(true_positive_rate + (args.nb_dec - decoy))/ float(args.nb_dec + args.nb_act)
        AUC = np.trapz(true_positive_rate_for_every_mol,false_positive_rate_for_every_mol)
        if args.verbose == True:
            print('******************')
            print('> Hits: ' + str( active + decoy))
            print('> True positive rate : ') + str(true_positive_rate)
            print('> False positive rate: ') + str(false_positive_rate)
            print('> Precision: ' + str(precision))
            print('> Accuracy: ' + str(accuracy))
            print('> Number of true positive hits: ' + str((active)))
            print('> Number of false positive hits: ' + str(decoy))
            print('> Accumulated area under the curve: ' + str(AUC))
            print('******************')

        return true_positive_rate_for_every_mol,false_positive_rate_for_every_mol


def draw_roc(true_positive_rate, false_positive_rate, args, auc_results):
    path,filename = os.path.split(str(args.input))
    fig = plt.figure(figsize = (6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(false_positive_rate,true_positive_rate)
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls='dotted', c='black')
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_xlabel('False Positive')
    ax.set_ylabel('True Positive')
    ax.set_title('ROC curve for ' + filename)

    ax.text(0.4,0.1,str(auc_results), bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
    fig.savefig(args.output)
    fig.clear()

def calculate_auc(true_positive_rate_for_every_mol, false_positive_rate_for_every_mol, args):
    TP_extrapolating = []
    FP_extrapolating = []
    treshold_list = [0.01, 0.02, 0.10, 0.20, 0.5, 1.0]
    nr_of_mol = args.nb_act + args.nb_dec

    for i in np.linspace(true_positive_rate_for_every_mol[-2], true_positive_rate_for_every_mol[-1], num = int(nr_of_mol) - len(true_positive_rate_for_every_mol), endpoint=True):
        TP_extrapolating.append(i)

    for i in np.linspace(false_positive_rate_for_every_mol[-2], false_positive_rate_for_every_mol[-1], num= int(nr_of_mol) -len(false_positive_rate_for_every_mol), endpoint=True):
        FP_extrapolating.append(i)

    TP_rate = true_positive_rate_for_every_mol[0:-2] + TP_extrapolating
    FP_rate = false_positive_rate_for_every_mol[0:-2] + FP_extrapolating
    auc_results = []
    for treshold in treshold_list:
        mol_treshold = float(nr_of_mol) * treshold

        for i in range(0, nr_of_mol):
            if FP_rate[i] >= treshold or TP_rate[i] >= treshold:
                boarder = max(FP_rate[i], TP_rate[i])
                FP_tmp = FP_rate[0:i+1] + [boarder]
                TP_tmp = TP_rate[0:i+1] + [boarder]
                area = np.trapz(TP_tmp, FP_tmp, dx=0.000001)
                scale_area = np.trapz([0,boarder], [0,boarder])*2
                result = np.round(area/scale_area, 2)
                auc_results.append(result)
                break

    if args.verbose == True:
        print( 'Calculating AUC values for :')
	legend = [treshold_list]
	legend.append('Nr of actives')
	legend.append('Nr of decoy')
   	print( legend)

    auc_results = [(np.round(float(x), 2)) for x in auc_results]
    auc_results = [(float(x)) for x in auc_results]
    auc_results.append(args.nb_act)
    auc_results.append(args.nb_dec)
    print( str(auc_results).replace('[', '').replace(']', ''))
    return auc_results

if __name__ == '__main__':

    args = parse_arguments()
    molecules = reading_input(args)
    true_positive_rate, false_positive_rate = calculating_TP_and_FP(molecules, args)

    auc_results = calculate_auc(true_positive_rate, false_positive_rate, args)
    if args.output:
        draw_roc(true_positive_rate, false_positive_rate, args, auc_results)
