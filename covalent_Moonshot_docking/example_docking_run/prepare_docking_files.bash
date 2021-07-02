#!/bin/bash

for f in *prot.pdb
	do prepare_receptor4.py -r $f -o ''${f::-4}'.pdbqt' -U 'nphs'	
	echo ''$f' done first step'
	done
	
for d in *flex.pdb
	do prepare_ligand4.py -l $d -R 1
	echo ''$d' done second step'
	done

for f in *_flex.pdbqt
	do cat flex_header.pdbqt $f flex_tail.pdbqt > ''${f::-6}'_clean.pdbqt'
	echo 'replaced header and tail'
	done
	
for f in *clean.pdbqt
	do prepare_dpf4.py -i template.dpf -l HOH.pdbqt -r  'Mpro-'${f:(-22):5}'_0_dim.sup_prep_prot.pdbqt' -x $f -o ${f::-16}'dock.dpf'
	echo ''$f' done dpf with Mpro-'${f:(-21):5}'_0_dim.sup_prep_prot.pdbqt' 
	done
	
for f in *clean.pdbqt	
	do prepare_gpf4.py -i template.gpf -l HOH.pdbqt -x $f -r 'Mpro-'${f:(-22):5}'_0_dim.sup_prep_prot.pdbqt' -o ${f::-16}'dock.gpf'
	echo ''$f' done gpf with Mpro-'${f:(-21):5}'_0_dim.sup_prep_prot.pdbqt' 
	done
	
for f in *dock.dpf
	do sed -i 's/unbound_model extended/unbound_model bound/' $f
	done
	
find . -name '*dock.gpf' | parallel --eta --joblog my_parallel.joblog /home/marc/custom_software/x86_64Linux3/autogrid4 -p {} -l {.}_log.glg 

find . -name '*dock.dpf' | parallel --eta --joblog my_parallel.joblog /home/marc/custom_software/x86_64Linux3/autodock4 -p {} -l {.}_log.dlg

