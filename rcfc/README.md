# RCFC

For an efficient implementation of the RCFC method, we want to stop the (primal) simplex algorithm as soon as the considered basic variable exceeds a value.
As a workaround, we called the simplex algorithm with a limited number of iterations and checked in between each run.
For that, we need to add a function to SCIP that yields the primal solution of a column, even when the simplex is not terminated yet.

The following steps need to be followed to use our code:
1. 
