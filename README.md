# Current datasets tested in PRS analysis:

PGCnoCLOZUK
-------------- 
Location on rocks: /home/SHARED/PGC/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt
Build: 37.13 : tested on one SNP, rs62513865 had the same BP position.
Summary stats

CLOZUK_no_PGC:
---------------
Location on rocks: /mnt/databank/CLOZUK/GWAS/BGE/*CLOZUK_GWAS_BGE*.tar.gz
Build: 
Best guess genotype data


# Preparing the config file
While not ideal, all preparations for the pipeline should be completed within _PRS\_arguments.sh_ and _PRS\_Pipeline\_for\_local.sh_

in _PRS\_arguments.sh_ where the line:
> elif [ "$whereami" = 'v1711-0ab8c3db.mobile.cf.ac.uk' ]; then

is found, replace **v1711-0ab8c3db.mobile.cf.ac.uk** with the character string found after running following command in a terminal:
> uname -n

for the home_OS argument, add the physical location of your hard-drive where the the analysis will be run to replace the root "**/**".
For example in MAC, most users have '/Users' while Unix/windows usually have '/home'. if unsure, open a terminal session and type:
> pwd

and extract the first directory after the root '**/**'

## Arguments
