/**
 * @file main.c
 * @brief Programme principal pour la simulation par éléments finis d'un rail UIC60.
 * @authors Edouard Meurant & Théodore Moulaert
 *
 * Ce programme effectue les étapes suivantes :
 * 1. Initialisation de la bibliothèque FEM et Gmsh.
 * 2. Génération, importation et correction du maillage.
 * 3. Attribution de noms aux frontières (domaines).
 * 4. Configuration du problème d'élasticité linéaire plane.
 * 5. Application des conditions aux limites (Dirichlet et Neumann).
 * 6. Résolution du système linéaire par éléments finis.
 * 7. Post-traitement des résultats (calculs, préparation visualisation).
 * 8. Visualisation interactive des résultats (déplacements, résidus) via OpenGL.
 * 9. Nettoyage et libération de la mémoire.
 */

 #include "glfem.h"
 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
 #include <stdbool.h>
 
 // =========================================================================
 //                            Helper Function
 // =========================================================================
 
 /**
  * @brief Fonction constante pour l'intégration de l'aire.
  */
 double funIntegrateArea(double x, double y) {
     return 1.0;
 }
 
 // =========================================================================
 //                              Main Program
 // =========================================================================
 
 int main(void) {
 
     // --- 0. User Instructions & Welcome Message ---
     printf("\n\n");
     printf("    FEM Simulation - Rail UIC60\n");
     printf("    ---------------------------\n");
     printf("    OpenGL Visualization Keys:\n");
     printf("      V : Displacement norm field (on deformed mesh)\n");
     printf("      D : Highlight domains (boundaries)\n");
     printf("      X : Horizontal residual forces field (on deformed mesh)\n");
     printf("      Y : Vertical residual forces field (on deformed mesh)\n");
     printf("      N : Cycle through highlighted domains\n");
     printf("      ESC : Exit visualization window\n");
     printf("\n\n");




 
     // =========================================================================
     //                      STAGE 1: Initialization
     // =========================================================================
     printf("STAGE 1: Initializing Geometry and FEM Library...\n");
 
     geoInitialize();
     femGeo* theGeometry = geoGetGeometry();
 
     theGeometry->h = 7.0;
     theGeometry->elementType = FEM_TRIANGLE;




 
     // =========================================================================
     //                 STAGE 2: Mesh Generation and Processing
     // =========================================================================
     printf("STAGE 2: Generating, Importing, and Fixing Mesh...\n");
 
     // Génération Gmsh, Importation dans femGeo, Correction en mémoire
     printf("  Generating mesh with Gmsh...\n");
     geoMeshGenerate();
     printf("  Importing mesh from Gmsh internal state...\n");
     geoMeshImport();
     printf("  Fixing mesh (removing unused nodes)...\n");
     geoMeshFix(theGeometry);
 
     // Vérification post-correction
     if (theGeometry->theNodes == NULL || theGeometry->theElements == NULL || theGeometry->theNodes->nNodes == 0) {
          fprintf(stderr,"Error: Mesh data is missing or empty after import/fix step.\n");
          geoFinalize();
          exit(EXIT_FAILURE);
     }
     printf("  Mesh imported and fixed. Nodes: %d, Elements: %d, Edges: %d, Domains: %d\n",
             theGeometry->theNodes->nNodes,
             theGeometry->theElements->nElem,
             (theGeometry->theEdges ? theGeometry->theEdges->nElem : 0),
             theGeometry->nDomains);




 
     // =========================================================================
     //                        STAGE 3: Domain Naming
     // =========================================================================
     printf("STAGE 3: Naming Domains (Boundaries)...\n");
 
     // !! ========================== ATTENTION CRUCIALE ========================== !!
     // !! Les indices iDomain (0, 1, 2, ...) dans geoSetDomainName DOIVENT       !!
     // !! correspondre à l'ORDRE dans lequel geoMeshImport a traité les          !!
     // !! Physical Groups de dimension 1 (lignes) de votre fichier .geo.         !!
     // !! Cet ordre N'EST PAS forcément le Tag (numéro) du Physical Group.       !!
     // !! MÉTHODE : Décommentez geoMeshPrint() ci-dessous pour vérifier l'ordre. !!
     // !! ====================================================================== !!
     // geoMeshPrint(); // Décommentez temporairement pour vérifier les indices/noms
 
     // Application des noms (VÉRIFIEZ ET CORRIGEZ LES INDICES !)
     if (theGeometry->nDomains > 0)  geoSetDomainName(0, "Bottom4");   else printf("  Warning: Skipping domain index 0 (nDomains=%d)\n", theGeometry->nDomains);
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
     if (theGeometry->nDomains > 22) geoSetDomainName(22, "Bottom0"); else if (theGeometry->nDomains <= 22) printf("  Warning: Skipping domain index 22\n");
 
     printf("  Domain naming attempt completed.\n");




 
     // =========================================================================
     //               STAGE 4: Optional Mesh Saving / Visu Setup
     // =========================================================================
     printf("STAGE 4: Optional Saving and Mesh Visualization Setup...\n");
 
     // Sauvegarde du maillage final pour inspection externe
     printf("  Writing final (fixed) mesh to elasticity_final.txt...\n");
     geoMeshWrite("../data/mesh.txt");
 
     // Calcul du champ de taille de référence pour visualisation
     femNodes *theNodes = theGeometry->theNodes;
     double *meshSizeField = malloc(theNodes->nNodes * sizeof(double));
     if (!meshSizeField) { Error("Failed to allocate meshSizeField"); }
     printf("  Calculating mesh size field for visualization...\n");
     for(int i=0; i < theNodes->nNodes; ++i) {
         meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
     }
 
     // Bloc de visualisation optionnelle du maillage
     printf("  Mesh visualization block present but disabled (enable by changing 'if(false)').\n");
     if (true) // Activer/déactiver la visualisation du maillage
     {
         printf("  Starting mesh visualization (optional)...\n");
         GLFWwindow* window_mesh = glfemInit("EPL1110 : Mesh generation (Fixed)");
         if (!window_mesh) { fprintf(stderr, "  Error: Failed to init GLFW window for mesh visualization\n"); }
         else {
             glfwMakeContextCurrent(window_mesh);
             int mode_mesh = 1;
             int domain_mesh = 0;
             int freezingButton_mesh = FALSE;
             double t_mesh, told_mesh = 0;
             char theMessage_mesh[256];
             do {
                 int w,h;
                 glfwGetFramebufferSize(window_mesh,&w,&h);
                 glfemReshapeWindows(theGeometry->theNodes,w,h);
                 t_mesh = glfwGetTime();
 
                 if (glfwGetKey(window_mesh, GLFW_KEY_D) == GLFW_PRESS) { mode_mesh = 0;}
                 if (glfwGetKey(window_mesh, GLFW_KEY_V) == GLFW_PRESS) { mode_mesh = 1;}
                 if (glfwGetKey(window_mesh, GLFW_KEY_N) == GLFW_PRESS && freezingButton_mesh == FALSE) {
                     if (theGeometry->nDomains > 0) domain_mesh = (domain_mesh + 1) % theGeometry->nDomains;
                     freezingButton_mesh = TRUE; told_mesh = t_mesh;
                 }
                 if (t_mesh-told_mesh > 0.5) {freezingButton_mesh = FALSE; }
 
                 glClearColor(1.0f,1.0f,1.0f,1.0f); glClear(GL_COLOR_BUFFER_BIT);
 
                 if (mode_mesh == 1) {
                     glfemPlotField(theGeometry->theElements, meshSizeField);
                     glfemPlotMesh(theGeometry->theElements);
                     sprintf(theMessage_mesh, "Mesh Size Field - Elements: %d ", theGeometry->theElements->nElem);
                 } else if (mode_mesh == 0) {
                     if (theGeometry->nDomains > 0 && domain_mesh < theGeometry->nDomains) {
                         if(theGeometry->theDomains[domain_mesh]) {
                              glfemPlotDomain( theGeometry->theDomains[domain_mesh]);
                              sprintf(theMessage_mesh, "Domain %d: %s (%d edges)", domain_mesh, theGeometry->theDomains[domain_mesh]->name, theGeometry->theDomains[domain_mesh]->nElem);
                         } else { sprintf(theMessage_mesh, "Domain %d: Invalid data", domain_mesh); }
                     } else { sprintf(theMessage_mesh, "No domains / Invalid index"); }
                 }
                  glColor3f(0.0f,0.0f,0.0f); glfemMessage(theMessage_mesh);
 
                 glfwSwapBuffers(window_mesh);
                 glfwPollEvents();
             } while( glfwGetKey(window_mesh,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
                     glfwWindowShouldClose(window_mesh) != 1 );
 
             glfwDestroyWindow(window_mesh);
             printf("  Mesh visualization finished.\n");
         }
     }




 
     // =========================================================================
     //                 STAGE 5: Elasticity Problem Setup
     // =========================================================================
     printf("STAGE 5: Setting up Elasticity Problem...\n");
 
     double E   = 211e9;
     double nu  = 0.3;
     double rho = 7850;
     double g   = 9.81;
 
     femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);
     if (!theProblem) {
         fprintf(stderr, "Error: femElasticityCreate failed!\n");
         free(meshSizeField);
         geoFinalize();
         exit(EXIT_FAILURE);
     }
     printf("  femProblem created. System size (DOFs): %d \n", theProblem->system->size);




 
     // =========================================================================
     //                 STAGE 6: Apply Boundary Conditions
     // =========================================================================
     printf("STAGE 6: Applying Boundary Conditions...\n");
 
     // Dirichlet X=0, Y=0 sur 'BottomX'
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
 
     // Dirichlet X=0, Y=0 sur 'GaucheX'
     femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_Y, 0.0);
 
     // Dirichlet X=0, Y=0 sur 'DroiteX'
     femElasticityAddBoundaryCondition(theProblem, "Droite0", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Droite0", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_Y, 0.0);
 
     // Dirichlet X=0, Y=0 sur 'HautX'
     femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut2", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut2", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_X, 0.0);
 
     // Neumann Fy = -1e4 N/m sur 'PatineX'
     femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, -1e4);
     femElasticityAddBoundaryCondition(theProblem, "Patine1", NEUMANN_Y, -1e4);
 
     // Affichage du résumé du problème après application des CL
     femElasticityPrint(theProblem);




 
     // =========================================================================
     //                 STAGE 7: Solve the FEM Problem
     // =========================================================================
     printf("STAGE 7: Assembling and Solving the Linear System...\n");
 
     double *theSoluce = femElasticitySolve(theProblem);
     printf("  System solved.\n");
     printf("  Calculating residual forces...\n");
     double *theForces = femElasticityForces(theProblem);
     printf("  Residual forces calculated.\n");
     printf("  Calculating area via integration...\n");
     double area = femElasticityIntegrate(theProblem, funIntegrateArea);
     printf("  Integration completed (Area = %e m^2).\n", area);




 
     // =========================================================================
     //         STAGE 8: Post-Processing and Visualization Setup
     // =========================================================================
     printf("STAGE 8: Post-Processing Results for Visualization...\n");
 
     double deformationFactor = 30000;
 
     // Allocation mémoire pour les champs nodaux de visualisation
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
     double *originalX = malloc(theNodes->nNodes * sizeof(double));
     double *originalY = malloc(theNodes->nNodes * sizeof(double));
      if (!normDisplacement || !forcesX || !forcesY || !originalX || !originalY) {
         free(normDisplacement); free(forcesX); free(forcesY); free(originalX);
         Error("Failed to allocate visualization/backup fields");
      }
 
     // Calcul des champs et application du facteur de déformation pour la visu
     printf("  Processing nodal results and applying deformation factor...\n");
     for (int i = 0; i < theNodes->nNodes; i++) {
         originalX[i] = theNodes->X[i];
         originalY[i] = theNodes->Y[i];
         if (theSoluce && theForces && (2*i + 1) < theProblem->system->size) {
             double ux = theSoluce[2 * i + 0];
             double uy = theSoluce[2 * i + 1];
             theNodes->X[i] += ux * deformationFactor;
             theNodes->Y[i] += uy * deformationFactor;
             normDisplacement[i] = sqrt(ux*ux + uy*uy);
             forcesX[i] = theForces[2 * i + 0];
             forcesY[i] = theForces[2 * i + 1];
         } else {
              fprintf(stderr, "  Warning: Error accessing solution/forces for node %d.\n", i);
              normDisplacement[i] = 0.0; forcesX[i] = 0.0; forcesY[i] = 0.0;
              theNodes->X[i] = originalX[i]; theNodes->Y[i] = originalY[i];
         }
     }




 
     // =========================================================================
     //                 STAGE 9: Output Numerical Results
     // =========================================================================
     printf("STAGE 9: Outputting Numerical Results...\n");
     geoMeshWrite("../data/elasticity.txt");
 
     double hMin = femMin(normDisplacement, theNodes->nNodes);
     double hMax = femMax(normDisplacement, theNodes->nNodes);
     printf("  ==== Minimum displacement norm     : %14.7e [m] \n", hMin);
     printf("  ==== Maximum displacement norm     : %14.7e [m] \n", hMax);
 
     double theGlobalForce[2] = {0.0, 0.0};
     if (theForces) {
         for (int i = 0; i < theNodes->nNodes; i++) {
             if ((2*i + 1) < theProblem->system->size) {
                 theGlobalForce[0] += theForces[2 * i + 0];
                 theGlobalForce[1] += theForces[2 * i + 1];
             }
         }
     }
     printf("  ==== Sum of Horizontal Residuals   : %14.7e [N] \n", theGlobalForce[0]);
     printf("  ==== Sum of Vertical Residuals     : %14.7e [N] \n", theGlobalForce[1]);
     printf("  ==== Estimated Weight (rho*g*Area*1m): %14.7e [N] \n", area * rho * g);




 
     // =========================================================================
     //                 STAGE 10: Results Visualization
     // =========================================================================
     printf("STAGE 10: Starting Results Visualization (OpenGL)...\n");
 
     GLFWwindow* window = glfemInit("EPL1110 : Linear Elasticity Results");
     if (!window) {
         fprintf(stderr, "  Error: Failed to initialize GLFW window for results visualization.\n");
     } else {
         glfwMakeContextCurrent(window);
         int mode = 1;
         int domain = 0;
         int freezingButton = FALSE;
         double t, told = 0.0;
         char theMessage[256];
 
         do {
             int w, h;
             glfwGetFramebufferSize(window, &w, &h);
             glfemReshapeWindows(theGeometry->theNodes, w, h);
             t = glfwGetTime();
 
             if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { mode = 0; }
             if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) { mode = 1; }
             if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) { mode = 2; }
             if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) { mode = 3; }
             if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS && freezingButton == FALSE) {
                 if (theGeometry->nDomains > 0) domain = (domain + 1) % theGeometry->nDomains;
                 freezingButton = TRUE; told = t;
             }
             if (t - told > 0.5) { freezingButton = FALSE; }
 
             glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
             glClear(GL_COLOR_BUFFER_BIT);
 
             switch(mode) {
                 case 0:
                     if (theGeometry->nDomains > 0 && domain < theGeometry->nDomains && theGeometry->theDomains[domain]) {
                         glfemPlotDomain(theGeometry->theDomains[domain]);
                         sprintf(theMessage, "Domain %d: %s (%d edges)", domain, theGeometry->theDomains[domain]->name, theGeometry->theDomains[domain]->nElem);
                     } else { sprintf(theMessage, "No domains or invalid index %d", domain); }
                     break;
                 case 1:
                     glfemPlotField(theGeometry->theElements, normDisplacement);
                     glfemPlotMesh(theGeometry->theElements);
                     sprintf(theMessage, "Displacement Norm (Min: %.2e, Max: %.2e) (Deformed Plot x%.0f)", hMin, hMax, deformationFactor);
                     break;
                 case 2:
                     glfemPlotField(theGeometry->theElements, forcesX);
                     glfemPlotMesh(theGeometry->theElements);
                     sprintf(theMessage, "Residual Forces X (Deformed Plot x%.0f)", deformationFactor);
                     break;
                 case 3:
                     glfemPlotField(theGeometry->theElements, forcesY);
                     glfemPlotMesh(theGeometry->theElements);
                     sprintf(theMessage, "Residual Forces Y (Deformed Plot x%.0f)", deformationFactor);
                     break;
                 default:
                      sprintf(theMessage, "Unknown display mode: %d", mode);
                      break;
             }
 
             glColor3f(0.0f, 0.0f, 0.0f);
             glfemMessage(theMessage);
 
             glfwSwapBuffers(window);
             glfwPollEvents();
 
         } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
                  glfwWindowShouldClose(window) != 1);
 
         glfwDestroyWindow(window);
         printf("  Results visualization finished.\n");
     }




 
     // =========================================================================
     //                 STAGE 11: Cleanup and Finalization
     // =========================================================================
     printf("STAGE 11: Cleaning Up Memory and Libraries...\n");
 
     // Restaurer les coordonnées nodales originales
     printf("  Restoring original node coordinates...\n");
     if (originalX && originalY) {
         for(int i=0; i < theNodes->nNodes; ++i) {
             theNodes->X[i] = originalX[i];
             theNodes->Y[i] = originalY[i];
         }
         free(originalX); originalX = NULL;
         free(originalY); originalY = NULL;
     }
 
     // Libérer la mémoire allouée dans main()
     printf("  Freeing allocated memory...\n");
     free(meshSizeField);    meshSizeField = NULL;
     free(normDisplacement); normDisplacement = NULL;
     free(forcesX);          forcesX = NULL;
     free(forcesY);          forcesY = NULL;
 
     // Libérer les structures du problème EF
     femElasticityFree(theProblem); theProblem = NULL;
 
     // Libérer la géométrie/maillage et finaliser Gmsh
     geoFinalize();
 
     // Terminer GLFW
     glfwTerminate();
 
     printf("Program finished successfully.\n");
     return 0;
 }