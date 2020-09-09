
# Molecular Dynamics runs

The molecular dynamics simulations can be started by making the md_run 
script executable and running it in the according folders:

```
chmod u+x md_run.sh
./md_run.sh
```

This script generates an folder output containin all of the run data of the 
MD. As described in the script this simulations contains of the usual steps of
energy minimization followed by equilibration in NVT and NPT ensemble, followed
by an unconstrained production run. The according configurations .MDP files
can be found in the input folder. Changes to the details of the according 
simulation steps like the number of MD steps or which algorithm is to be used
for minimization, need to be done by changing these .MDP files.

For using the script one may need to adapt the first lining sourcing gromacs 

```
source /opt/gromacs-2019/bin/GMXRC
```

This needs to be adapted to path where gromacs is located. In general it is
also sufficient to just have gmx or gmx_mpi somewhere in the PATH variable.
Meaning we could also do something like this 

```
export PATH=/opt/gromacs-2019/bin/GMXRC:$PATH
```


