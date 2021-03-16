# runs iqtree to generate ancestral state of each node per patient.


nohup iqtree -s data/tree/15664.fasta -m GTR+G4+FO -bb 1000 -alrt 1000 -asr -nt 6

nohup iqtree -s data/tree/16207.fasta -m GTR+G4+FO -bb 1000 -alrt 1000 -asr -nt 6

nohup iqtree -s data/tree/22763.fasta -m GTR+G4+FO -bb 1000 -alrt 1000 -asr -nt 6