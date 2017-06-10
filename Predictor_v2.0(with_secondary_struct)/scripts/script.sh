#!/bin/bash 
 grep "^SEQUENCE" "../prediction/3V7E_C_5FK3_A/files/emboss.aln" | sed -e "s/^SEQUENCE/3V7E_C/;N" | sed -e "s/^SEQUENCE/5FK3_A/" | sed -e "s/ [[:digit:]]*/ /g" > ../prediction/3V7E_C_5FK3_A/files/tta.aln
grep "^ATOM.................A" ../pdbs/5fk3.pdb > ../prediction/3V7E_C_5FK3_A/files/template.pdb