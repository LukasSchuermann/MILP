# RCFC

For an efficient implementation of the RCFC method, we want to stop the (primal) simplex algorithm as soon as the considered basic variable exceeds a value.
As a workaround, we called the simplex algorithm with a limited number of iterations and checked in between each run.
For that, we need to add a function to SCIP that yields the primal solution of a column, even when the simplex is not terminated yet.

The following steps need to be followed to use our code:
  1. Download the source code of the SCIP Optimization Suite from https://scipopt.org/index.php#download (we used version 9.1.1)
  2. In "scip/src/" we need to adapt "scip/lpi.h" and the corresponding LP-solver file ("lpi/lpi_spx2" in our case) accordingly:
     ```markdown
     # Add this function to scip/lpi.h
     SCIP_EXPORT
     SCIP_Real SCIPlpiColGetNewLPvalRCFC(
        SCIP_LPI*           lpi,
        int                 colIndex
     );
     # Add this declaration to lpi/lpi_spx2
     SCIP_Real SCIPlpiColGetNewLPvalRCFC(
        SCIP_LPI*           lpi,
        int                 colIndex
     ){
        return lpi->spx->getPrimalRealIndex(colIndex);
     }
     ```
  3. In "soplex/src/soplex/" we need to add the following code to the corresponding files:
     ```markdown
     # Add this function to the SoPlexBase class in soplex.h (e.g. at line 664)
     double getPrimalRealIndex(int index);
     # Add the declaration to soplex.hpp
     template <class R>
     double SoPlexBase<R>::getPrimalRealIndex(int index){
        return _solReal._primal[index];
     }
     # Add this function to soplex_interface.h
     double SoPlex_getPrimalRealIndex(void* soplex, int index);
     # Add the declaration to soplex_interface.cpp
     double SoPlex_getPrimalRealIndex(void* soplex, int index){
        SoPlex* so = (SoPlex*)(soplex);
        return so->getPrimalRealIndex(index);
     }
     ```
