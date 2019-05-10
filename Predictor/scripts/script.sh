#!/bin/sh 
 grep "^SEQUENCE" "../prediction/4TNA_A_486D_A/files/emboss.aln" | sed -e "s/^SEQUENCE/4TNA_A/;N" | sed -e "s/^SEQUENCE/486D_A/" | sed -e "s/ [[:digit:]]*/ /g" > ../prediction/4TNA_A_486D_A/files/tta.aln
grep "^ATOM.................A" ../pdbs/486d.pdb > ../prediction/4TNA_A_486D_A/files/template.pdb