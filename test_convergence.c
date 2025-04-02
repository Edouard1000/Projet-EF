#include "glfem.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double TARGET_POINT_X_29 = 19.4;
const double TARGET_POINT_Y_29 = 171.0;


int main(void) {
    FILE *file = fopen("Gaussian_elimination_convergence.csv", "w");
    if (!file) { perror("fopen"); exit(EXIT_FAILURE); }
    fprintf(file, "h,solve_time\n");

    for (double h = 30.0; h >= 20.0; h -= 0.5) {
     printf("BOUCLE SUR h = %f\n", h);
              
 
     // =========================================================================
     //                      STAGE 1: Initialization
     // =========================================================================
     printf("STAGE 1: Initializing Geometry and FEM Library...\n");
     geoInitialize();
     femGeo* theGeometry = geoGetGeometry();
     theGeometry->h = h; // Valeur de h pour cette simulation unique
     theGeometry->elementType = FEM_TRIANGLE;
 
 
     // =========================================================================
     //                 STAGE 2: Mesh Generation and Processing
     // =========================================================================
     printf("STAGE 2: Generating, Importing, and Fixing Mesh...\n");
     geoMeshGenerate();
     geoMeshImport();
     geoMeshFix(theGeometry);
     if (!theGeometry->theNodes || !theGeometry->theElements || theGeometry->theNodes->nNodes == 0) {
          fprintf(stderr,"Error: Mesh data is missing or empty after import/fix step.\n");
          geoFinalize(); exit(EXIT_FAILURE);
     }
     printf("  Mesh imported and fixed. Nodes: %d, Elements: %d, Edges: %d, Domains: %d\n",
             theGeometry->theNodes->nNodes, theGeometry->theElements->nElem,
             (theGeometry->theEdges ? theGeometry->theEdges->nElem : 0), theGeometry->nDomains);
 
 
     // =========================================================================
     //                        STAGE 3: Domain Naming
     // =========================================================================
     printf("STAGE 3: Naming Domains (Boundaries)...\n");
     // Votre bloc de geoSetDomainName original (copiez-le ici)
     if (theGeometry->nDomains > 0)  geoSetDomainName(0, "Bottom4");   else printf("  Warning: Skipping domain index 0\n");
     if (theGeometry->nDomains > 1)  geoSetDomainName(1, "Bottom5");
     if (theGeometry->nDomains > 2)  geoSetDomainName(2, "Bottom6");
     if (theGeometry->nDomains > 3)  geoSetDomainName(3, "Bottom7");
     if (theGeometry->nDomains > 4)  geoSetDomainName(4, "Bottom1");
     if (theGeometry->nDomains > 5)  geoSetDomainName(5, "Bottom2");
     if (theGeometry->nDomains > 6)  geoSetDomainName(6, "Bottom3");
     if (theGeometry->nDomains > 7)  geoSetDomainName(7, "Droite0");
     if (theGeometry->nDomains > 8)  geoSetDomainName(8, "Gauche0");
     if (theGeometry->nDomains > 9)  geoSetDomainName(9, "Droite1");
     if (theGeometry->nDomains > 10) geoSetDomainName(10, "Gauche1");
     if (theGeometry->nDomains > 11) geoSetDomainName(11, "Droite2");
     if (theGeometry->nDomains > 12) geoSetDomainName(12, "Haut0");
     if (theGeometry->nDomains > 13) geoSetDomainName(13, "Haut5");
     if (theGeometry->nDomains > 14) geoSetDomainName(14, "Gauche2");
     if (theGeometry->nDomains > 15) geoSetDomainName(15, "Patine0");
     if (theGeometry->nDomains > 16) geoSetDomainName(16, "Patine1");
     if (theGeometry->nDomains > 17) geoSetDomainName(17, "Bottom8");
     if (theGeometry->nDomains > 18) geoSetDomainName(18, "Haut4");
     if (theGeometry->nDomains > 19) geoSetDomainName(19, "Haut3");
     if (theGeometry->nDomains > 20) geoSetDomainName(20, "Haut2");
     if (theGeometry->nDomains > 21) geoSetDomainName(21, "Haut1");
     if (theGeometry->nDomains > 22) geoSetDomainName(22, "Bottom0"); else printf("  Warning: Skipping domain index 22\n");
     printf("  Domain naming attempt completed.\n");
 
 
     
      // =========================================================================
      //                 STAGE 5: Elasticity Problem Setup
      // =========================================================================
      printf("STAGE 5: Setting up Elasticity Problem...\n");
  
      double E   = 211e9;
      double nu  = 0.3;
      double rho = 10000.0;
      double g   = 9.81;
  
      femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);
      if (theProblem == NULL || theProblem->system == NULL) {
        fprintf(stderr, "Erreur : theProblem ou system est NULL !\n");
        geoFinalize(); // toujours fermer proprement
        exit(EXIT_FAILURE);
     }
      printf("  femProblem created. System size (DOFs): %d \n", theProblem->system->size);
 
 
 
     // =========================================================================
     //                 STAGE 6: Apply Boundary Conditions
     // =========================================================================
     printf("STAGE 6: Applying Boundary Conditions...\n");
     
    
    femElasticityAddBoundaryCondition(theProblem, "Bottom0", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom2", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom3", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom4", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom5", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom6", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom7", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom8", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom0", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom2", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom3", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom4", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom5", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom6", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom7", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom8", DIRICHLET_X, 0.0);
    // Gauche
    //femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_X, 0.0);
    //femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);
    //femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_Y, 0.0);
    //femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_Y, 0.0);
    // Droite
    //femElasticityAddBoundaryCondition(theProblem, "Droite0", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_X, 0.0);
    //femElasticityAddBoundaryCondition(theProblem, "Droite0", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_Y, 0.0);
    // Haut
    //femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_Y, 0.0); 
    femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_X, 0.0);
    //femElasticityAddBoundaryCondition(theProblem, "Haut2", DIRICHLET_Y, 0.0); 
    // Charge Neumann
    double load = -6250000; // GPa

    femElasticityAddBoundaryCondition(theProblem, "Haut2", NEUMANN_Y, load*0.5);
    femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, load*0.5);
    femElasticityAddBoundaryCondition(theProblem, "Haut1", NEUMANN_Y, load*0.5);
        

    femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_Y, 0.0);
     
     femElasticityPrint(theProblem);
 
 
     // =========================================================================
     //                 STAGE 7: Solve the FEM Problem
     // =========================================================================
     printf("STAGE 7: Assembling and Solving the Linear System...\n");
     
     clock_t start = clock();
     double *theSoluce = femElasticitySolve(theProblem);
     clock_t end = clock();
     printf("  System solved.\n");

     double solveTime = (double)(end - start) / CLOCKS_PER_SEC;

     printf("  Calculating residual forces...\n");
     double *theForces = femElasticityForces(theProblem);

     
       

     fprintf(file, "%.2f,%.6f\n", h, solveTime);

     femElasticityFree(theProblem); 
     geoFinalize();
    }

    fclose(file);
    printf("\n>>> Résultats sauvegardés dans convergence_results.csv \n");
    return 0;
}
