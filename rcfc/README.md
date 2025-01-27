# RCFC

For an efficient implementation of the RCFC method, we want to stop the (primal) simplex algorithm as soon as the considered basic variable exceeds a value.
As a workaround, we called the simplex algorithm with a limited number of iterations and checked in between each run.
For that, we need to add a function to SCIP that yields the primal solution of a column, even when the simplex is not terminated yet.

The following steps need to be followed to use our code:
  1. Download the source code of the SCIP Optimization Suite from https://scipopt.org/index.php#download (we used version 9.1.1)
  2. In "scip/src/" we need to adapt "scip/lpi.h" and the corresponding LP-solver file ("lpi/lpi_spx2" in our case) accordingly:
     ```markdown
     SCIP_EXPORT
     SCIP_Real SCIPluklpiColGetNewLPval(
     SCIP_LPI*           lpi,
     int                 colIndex
     );
     ```
