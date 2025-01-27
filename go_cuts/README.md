# GO cuts

This project contains the implementation of the GO cuts separator for MILPs of my PhD thesis.
It is implementend in C++ using the MINLP solver SCIP.

## Requirements

To run, this codes needs an installation of SCIP (we used version 9.1.1). Env_Variable "SCIP_DIR" shall contain the path to the SCIP installation.
See the website for an installation guide:
https://scipopt.org/#scipoptsuite

```markdown
$ mkdir build
$ cd build
$ cmake .. -DSCIP_DIR="/path/to/scip/installation"
$ make
```
Now we can run the code as in the following example:
```markdown
$ ./go ../../instances/pure_integer/irp.mps +useCutting
```
