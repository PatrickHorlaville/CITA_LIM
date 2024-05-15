This code reads in a halo catalogue, applies a HOD (Shang or Manera) and 
projects a CIB map with the Planck 2013 or 2015 model or an optical map with 
the Manera et al. 2013 model

To get code: git clone git@gitlab.com:georgestein/hod2map.git

Main code driver is hod2map.py. Just use mpirun -np <nproc> ./hod2map.py 
Maps are saved to maps/

The massive number of satellites causes a large memory footprint. 
For this reason large jobs for now need to be submitted to a 32Gb node. 
A sample qsub script can be seen in submit_map.sh. 
ie include the line:
#PBS -l nodes=1:m32g:ppn=8,walltime=1:00:00


Parameter file - params.py
All parameters you would need to change are in here, 
except for the observation frequencies which are hardcoded in the driver.


Modules needed to run hod2map are:
module load gcc intel/13.1.1 python/2.7.5
module load openmpi/intel/1.6.4
