#! /bin/ksh
# @ job_name   = atom_dz
# @ job_type=MPICH
# @ error = job.err.$(jobid)
# @ output = job.out.$(jobid)
# @ class = n002
# @ restart = no
# @ node = 2
# @ tasks_per_node = 20
# @ wall_clock_limit = 00:10:00
# @ notification = complete
# @ notify_user = sharma@fkf.mpg.de
# @ queue

#. /usr/share/modules/init/bash
#export TMPDIR=/ptmp/avl/

module load ifort/14.0.2
module load mpi.intel/4.1.3

mpiexec.hydra ./block.spin_adapted dmrg.conf >dmrg.out
#mpiexec.hydra ./block.spin_adapted compress.conf >compress.out
#mpiexec.hydra ./block.spin_adapted response.conf >response.out
#mpiexec.hydra ./a.out >out

