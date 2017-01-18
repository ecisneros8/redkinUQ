#!/bin/bash
################################################################
##                                                            ##
##                    Campus Cluster                          ##
##           Sample MVAPICH2 Job Batch Script                 ##
##                                                            ##
## PBS Options                                                ##
##                                                            ##
##  option -l                                                 ##
##    walltime: maximum wall clock time (hh:mm:ss)            ##
##       nodes: number of 12-core(taub) nodes                 ##
##                        16-core(golub) nodes                ##
##                        20-core(golub) nodes                ##
##         ppn: cores per node to use (1 thru 12) on taub     ##
##                                    (1 thru 16) on golub    ##
##                                    (1 thru 20) on golub    ##
##                                                            ##
##  option -N                                                 ##
##    job name (default = name of script file)                ##
##                                                            ##
##  option -q                                                 ##
##    queue name ( -q name_of_queue )                         ##
##                                                            ##
##  option -o                                                 ##
##     filename for standard output at end of job             ##
##     (default = <job_name>.o<job_id>).  File is written to  ##
##     directory from which qsub was executed. Remove extra   ##
##     "##" from the PBS -o line if you want to name your     ##
##     own file.                                              ##
##                                                            ##
##  option -e                                                 ##
##     filename for standard error at end of job              ##
##     (default = <job_name>.e<job_id>).  File is written to  ##
##     directory from which qsub was executed. Remove extra   ##
##     "##" from the PBS -e line if you want to name your     ##
##     own file.                                              ##
##                                                            ##
##  option -j                                                 ##
##     Join the standard output and standard error streams    ##
##     of the job                                             ##
##     ( -j oe  merges stderr with stdout and generates a     ## 
##              stdout file <job_name>.o<job_id>              ##
##       -j eo  merges stdout with stderr and generates a     ##
##              stderr file <job_name>.e<job_id>  )           ##
##                                                            ##
##  option -m                                                 ##
##     mail_options (email notifications)                     ##
##     The mail_options argument is a string which consists   ## 
##     of either the single character "n", or one or more of  ##
##     the characters "a", "b", and "e".                      ##
##     ( -m a   Send mail when the job is aborted.            ##
##       -m be  Send mail when the job begins execution and   ##
##              when the job terminates.                      ##
##       -m n   Do not send mail.  )                          ##
##                                                            ##
################################################################
#
#PBS -l walltime=00:10:00
#PBS -l nodes=1:ppn=20
#PBS -N RCCE
#PBS -q cse 
#PBS -j oe
#PBS -o rcce.out                                                            
#PBS -e rcce.err
###PBS -m be
#
#####################################

# Change to the directory from which the batch job was submitted
cd $PBS_O_WORKDIR

# Run the MPI code
mpiexec ./redkin $NC $JOB

