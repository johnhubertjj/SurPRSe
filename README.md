# Super-PRS

* A workflow in Bash, R and python designed to allow easy processing and production of polygenic risk scores (both whole genome and set scores) on [ARCCA's](http://www.cardiff.ac.uk/advanced-research-computing) Raven supercomputer
* For a [tutorial](http://gitlab.psycm.cf.ac.uk/john/Schizophrenia_PRS_pipeline_scripts/wikis/example-prs-tutorial), [debugging](http://gitlab.psycm.cf.ac.uk/john/Schizophrenia_PRS_pipeline_scripts/wikis/debugging_long) information and other documentation, please refer to the [wiki](http://gitlab.psycm.cf.ac.uk/john/Schizophrenia_PRS_pipeline_scripts/wikis/home)
* If used for the creation of results to be published (one can dream...), please acknowledge  John Hubert <a href="mailto:hubertjj@cardiff.ac.uk">hubertjj@cardiff.ac.uk</a> and Valentina Escott-Price <a href="mailto:escottpricev@cardiff.ac.uk">escottpricev@cardiff.ac.uk</a>
* I'm always looking for people to help in making this better, If you are familiar to git please fork this repository. 
    * The code has a General Public license so you are free to make changes of your own just as long as you document them!
    * If you are unfamiliar with the 'git' structure but still want to help out, please email John Hubert @ <a href="mailto:hubertjj@cardiff.ac.uk">hubertjj@cardiff.ac.uk</a> and we can sort something out!
* If you notice a bug or a fix that is required please add to the [issues page](http://gitlab.psycm.cf.ac.uk/john/Schizophrenia_PRS_pipeline_scripts/issues) or email <a href="mailto:hubertjj@cardiff.ac.uk">hubertjj@cardiff.ac.uk</a>  
* Please note that 'Super' refers to Supercomputer and not the state of the workflow, "Okay-PRS" didn't really have the same 'kick' to it...  

 
# SAVE-PRS

** COMING SOON**

* [**S**uper-PRS **A**nalysis **V**iewing **E**nvironment - **P**olygenic **R**isk **S**cores](https://johnhubertjj.shinyapps.io/Viewing_PRS_two_files/)
* Comparing whole genome polygenic risk scores to polygenic risk set scores is difficult as the 'ideal' significance threshold for a whole genome polygenic risk score is unlikely to be the same for a polygenic risk set score, and will likely change depending on which set you test.
* SAVE-PRS allows for "real-time" comparison of set scores to whole genome polygenic risk scores in a [shiny app](https://shiny.rstudio.com/), to allow comparison across all significance thresholds and allow better interpretation of all results.

# Current datasets tested in PRS analysis:

PGC2 (SCZ)
-------------- 

#### PGC2minusCLOZUK (SCZ2minusCLOZUK)
Location on rocks: /home/SHARED/PGC/daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_noclo.txt  
Summary stats

#### Sweden+SCZ1
Location: https://www.med.unc.edu/pgc/results-and-downloads  
Summary stats  

#### SCZ1
Location: https://www.med.unc.edu/pgc/results-and-downloads  
Summary stats  
**WARNING** build is UCSC hg18, will likely need to convert to hg19.


CLOZUK
---------------

#### CLOZUK_no_PGC
Location on raven: /neurocluster/databank/CLOZUK/GWAS/BGE/\*CLOZUK_GWAS_BGE\*.tar.gz  
Best guess genotype data

#### CLOZUK_Meta-analysis with PGC
Location on raven: /neurocluster/databank/CLOZUK/GWAS/SUMSTATS/CLOZUK_PGC2noclo.wCHRX.w1000Gfrq.METAL.assoc.dosage.gz  
Summary stats

ALSPAC
---------------
Best guess genotype data

Biobank wv 1
---------------
Best guess genotype data

