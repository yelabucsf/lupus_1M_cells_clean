#!/bin/bash                        
#                           
#$ -S /bin/bash                  
#$ -o /netapp/home/mincheol/outputs
#$ -e /netapp/home/mincheol/outputs
#$ -cwd                          
#$ -r y                          
#$ -j y                          
#$ -l mem_free=20G             
#$ -l h_rt=3:00:00

source activate scvi
python /netapp/home/mincheol/lupus_paper/figure3.myeloid/inteferon_axis.py
source deactivate