 #!/bin/bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/sharcnet/acml/4.4.0/gfortran-64bit/gfortran64/lib
module unload mkl intel
module load acml/gfortran-int64/4.4.0
export PATH=$PATH:/work/kimt33/sources/peter_oo/
