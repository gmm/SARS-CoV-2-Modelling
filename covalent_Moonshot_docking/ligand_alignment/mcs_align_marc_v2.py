#!/usr/bin/env python

# inspiration see https://sourceforge.net/p/rdkit/mailman/message/30641401/
# Updated to Python 3 and RDKit 2020.03.1 by Garrett M. Morris and Marc Moesser
# Expanded by Marc Moesser

import sys

if (len(sys.argv) < 4):
    print('MCSalign.py [template SDfile] [library SDF] [output SDF] [rms-threshold=-1] [shape-threshold=2] [nconfs=100] [matrix_check=True]')
    exit(1)
    
    
#new imports by marc
import scipy.spatial as ss
import numpy as np
import pandas as pd


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolAlign

# template for alignment
master = Chem.SDMolSupplier(sys.argv[1])
# database of 2D (or 3D) ligands to be aligned
suppl = Chem.SDMolSupplier(sys.argv[2])
# output-file
writer = Chem.SDWriter(sys.argv[3])
# prewriter = Chem.SDWriter('preopt.sdf') # to be removed

refmol = master[0]  # currently only the first cpds is used as reference -> to be exanded
if (refmol is None):
    print(f'Error for reference-molecule {sys.argv[1]}')
    exit(13)

#--------
#RMS threshold and shapethreshold are set too low so it is never reachable. The script will therefore only
#output the best conformation for each molecule instead of all mols that pass the test. 
#Adjust these to your liking. Initial experiments point towards 0.8, 0.8 to be good values for both.
#matrix_check is a new functionality introduced by Marc Moesser to check for unrealistic conformations that have
#stacked atoms or unrealistic bonds. Only checks for too small distances and does not correct for large distances.
# - MARC MOESSER
#--------

rmsthreshold = -1
shapethreshold = -1
nconf = 100
matrix_check = True

if (len(sys.argv) > 4):
    rmsthreshold = float(sys.argv[4])
if (len(sys.argv) > 5):
    shapethreshold = float(sys.argv[5])
if (len(sys.argv) > 6):
    nconf = int(sys.argv[6])
if (len(sys.argv) > 7):
    matrix_check = int(sys.argv[7])

#
forcetol = 0.1
# rmsd for intra-comparisons and conf.gen. pruning cutoff
cutoff = 0.5

for i, mol in enumerate(suppl):
    if mol is None:
        print(f'Error for molecule at position {i:d}')  # skip bad molecules
    else:
        s = f'molecule {i:d}:'
        if (mol.HasProp('_Name')):
            s = f'{mol.GetProp("_Name")}'

        refconf = refmol.GetConformer()
        molconf = mol.GetConformer()
        mols = [refmol, mol]

        nsuccess = 0


        # NEW ADDITION: completeRingsOnly=True!!! this will ensure that generated conformers don't have forced,
        # impossible angles.    - MARC MOESSER

        mcsres = rdFMCS.FindMCS(mols, timeout=5, completeRingsOnly=True)
        patt = Chem.MolFromSmarts(mcsres.smartsString)

        #------------------------------
        # NEW ADDITION: Logic checks to see if the target molecule is somehow flawed and shouldnt be used in this script
        # - MARC MOESSER

        #Check for common substructure between the found mcs pattern and cysteine
        cysteine = 'N[C@H](C=O)CSC'
        cys = Chem.MolFromSmiles(cysteine)
        cys_test_mols = [cys, patt]
        mcscys = rdFMCS.FindMCS(cys_test_mols, timeout=5)

        #Reject mol if mcs its too small
        if mcsres.numAtoms < 3 or mcsres.canceled:
            print(f'Could not generate sufficient MCS for {s}')
        #Reject mol if mcs does not contain the cysteine
        elif mcscys.numAtoms < 6:
            print(f'Generated MCS of {s} does not contain Cys - ABORTED')
        #--------------------------------

        else:
            refat = refmol.GetSubstructMatch(patt)
            dbat = mol.GetSubstructMatch(patt)
            aliatoms = []
            coordmap = {}

            # generate atom-idx mapping  using coordinates from reference molecule
            for k in range(0, len(refat)):
                pt3 = refconf.GetAtomPosition(refat[k])
                molconf.SetAtomPosition(dbat[k], pt3)
                coordmap[dbat[k]] = pt3
                aliatoms.append([dbat[k], refat[k]])

            mol.AddConformer(molconf, assignId=0)
            vec = AllChem.EmbedMultipleConfs(mol, coordMap=coordmap, forceTol=forcetol, numConfs=nconf,
                                             pruneRmsThresh=cutoff)

            if (len(vec) == 0):
                print(f'No intial conformations generated for molecule {s}')
            else:
                bestsim = 1
                bestrms = 1e10
                conf = -1

                check = []
                # loop over all generated conformations
                for c in range(len(vec)):
                    check.append(1)
                    # optimize current conformation
                    ## @@ Original code used AllChem.UFFOptimizeMolecule(mol,confId=c) #,confId=0 , ,maxIters=100
                    # --- GMM updated code with this snippet: @@

                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=c)

                    # distance constraints for minimization
                    for k in range(len(refat)):
                        pt3 = refconf.GetAtomPosition(refat[k])
                        pIdx = ff.AddExtraPoint(pt3.x, pt3.y, pt3.z, fixed=True) - 1
                        ff.AddDistanceConstraint(pIdx, dbat[k], 0, 0, 100.)

                    #-----------------
                    #This is the crucial new line: ALIGN MOLECULES BEFORE ENERGY OPT!!! Otherwise some bonds
                    #may have impossible angles!!! - MARC MOESSER
                    rms_pre = Chem.rdMolAlign.AlignMol(mol, refmol, atomMap=aliatoms, prbCid=c)
                    #-----------------

                    ff.Initialize()
                    n = 4
                    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
                    while more and n:
                        more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
                        n -= 1
                    # --- GMM snippet ends --- @@

                    # compare current conformation with previously optimized conformations
                    # and ignore if it is to similar to any other conformation

                    for d in range(c):
                        if check[c] != 0 and c != d:
                            rms = AllChem.AlignMol(mol, mol, c, d)
                            if rms <= cutoff:
                                check[c] = 0
                                break

                    if check[c] == 1:
                        # re-align current conformation
                        rms = Chem.rdMolAlign.AlignMol(mol, refmol, atomMap=aliatoms, prbCid=c)
                        sim = AllChem.ShapeTanimotoDist(refmol, mol, confId2=c)

                        if matrix_check:
                        #CHECK FOR FUNKY ATOM STACKING ----- MARC MOESSER

                        # This block will calculate a distance matrix between all atoms in the conformer
                            atnum = mol.GetNumAtoms()
                            conf_to_check = mol.GetConformer(c)
                            pt_list = []
                            for j in range(0, atnum):
                                coord = conf_to_check.GetAtomPosition(j)
                                pt_list.append(coord)
                            #calc distance matrix in Angstrom between all atoms in conf
                            matrix = ss.distance_matrix(pt_list, pt_list)

                            #put in df, delete all 0 to search for values < than distance threshold. here 0.9 A.
                            df_temp = pd.DataFrame(matrix)
                            df_temp = df_temp.replace(0, np.nan)
                            if df_temp[(df_temp < 0.9).any()].empty:    
                                # THIS VALUE (0.9) IS ARBITRARY:
                                # <0.9 chosen since a normal bond should be at least 1 A.
                                # write out all conformations which are below both thresholds
                                # Change check to sim "<" shapethreshold since small sim == more similar - MARC MOESSER
                                if (rms < rmsthreshold and sim < shapethreshold):
                                    print(f'Writing conformation with RMS {rms:f} SIM {sim:f}')
                                    writer.write(mol, confId=c)
                                    nsuccess = nsuccess + 1

                                # if sim, rmsd are not below threshold -> store best solution so far
                                if sim < bestsim and rms < bestrms:
                                    bestsim = sim
                                    bestrms = rms
                                    conf = c
                            else:
                                print("rejected by distance matrix")
                        
                        else:
                            # write out all conformations which are below at least one threshold
                            #Change check to sim "<" shapethreshold since small sim == more similar - MARC MOESSER
                            if (rms < rmsthreshold or sim < shapethreshold):
                                print(f'Writing conformation with RMS {rms:f} SIM {sim:f}')
                                writer.write(mol, confId=c)
                                nsuccess = nsuccess + 1

                            # if sim, rmsd are not below threshold -> store best solution so far
                            if sim > bestsim and rms < bestrms:
                                bestsim = sim
                                bestrms = rms
                                conf = c

                if nsuccess == 0 and conf >= 0 and check[conf] != 0:
                    # only write best conformation if none has been saved so far
                    # This check looks for unrealistic thrsehold to output only the best compound. - MARC MOESSER
                    if rmsthreshold < 0 and shapethreshold < 0:
                        print(f'Writing best conformation shape-sim {bestsim:.3f} rms {bestrms:.3f}')
                        writer.write(mol, confId=conf)
                        nsuccess = 1
                if nsuccess == 0:
                    print('Error when generating conformation')

        print(f'Wrote {nsuccess:d} conformations for {s}')
