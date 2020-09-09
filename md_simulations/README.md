
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

## Evaluation

For evaluation of the according MD runs there are two scripts placed in each
of the MD run folders. First there ist remove_boundary, which is a bash
script using gmx trjconv to remove the periodic boundary conditions as well as 
water from the resulting trajectory. Also the structure is fitted with respect
to rotational and translational degrees of freedom, such that drifts and global 
rotations of the Molecule are corrected in the resulting trajectory.

Futhermore there is a python script evaluate_new_md.py or evaluate_md.py that
uses the biotite package for generating basic plots of RMSD, RMSF as well as
radius of gyration.

As a note: It may be necessary to change gmx to gmx_mpi in the remove_trajectoy 
scripts, or the other way around. However, this dependends on how
you gromacs version was compiled and if you want to use the mpi only version or
not.


