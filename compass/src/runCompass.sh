# Description: run Compass, requires a running Compass environment, see https://yoseflab.github.io/Compass/install.html for install instructions
#
# PRE: SCT matrix was generated with Seurat and log.data was reverted (see transform2delog.R)
#
# OUT: generates reactions.tsv in the respective output folders
#
# author Sascha Schäuble
# date of creation: Tue June 12 12:51:07 2024
# license: MIT

# brain
compass --data ../dat/snRNA/brain_SCT_matrix_delog.tsv --select-subsystems ../dat/modelDat/poi.txt --num-processes 60 --species mus_musculus --output-dir ../res/scores/br/
# kidney
compass --data ../dat/snRNA/kidney_SCT_matrix_delog.tsv --select-subsystems ../dat/modelDat/poi.txt --num-processes 60 --species mus_musculus --output-dir ../res/scores/ki/
# liver
compass --data ../dat/snRNA/liver_SCT_matrix_delog.tsv --select-subsystems ../dat/modelDat/poi.txt --num-processes 60 --species mus_musculus --output-dir ../res/scores/li/
# WAT
compass --data ../dat/snRNA/wat_SCT_matrix_delog.tsv --select-subsystems ../dat/modelDat/poi.txt --num-processes 60 --species mus_musculus --output-dir ../res/scores/wa/

