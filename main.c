/**
 * @file main.c
 * @brief Programme principal pour la simulation par éléments finis d'un rail UIC60.
 * @authors Edouard Meurant & Théodore Moulaert
 *
 * Ce programme effectue les étapes suivantes :
 * 1. Initialisation de la bibliothèque FEM et Gmsh.
 * 2. Génération, importation et correction du maillage.
 * 3. Attribution de noms aux frontières (domaines).
 * 4. Sauvegarde optionnelle et visualisation du maillage corrigé.
 * 5. Configuration du problème d'élasticité linéaire plane.
 * 6. Application des conditions aux limites (Dirichlet et Neumann).
 * 7. Résolution du système linéaire par éléments finis.
 * 8. Calcul des forces résiduelles, intégration de l'aire et récupération du déplacement du point cible.
 * 9. Post-traitement des résultats pour la visualisation (déplacements, forces).
 * 10. Visualisation interactive des résultats (déplacements, résidus) via OpenGL.
 * 11. Nettoyage et libération de la mémoire.
 */

 #include "glfem.h"
 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
 #include <stdbool.h>
 
 // =========================================================================
 //                            Target Point Definition
 // =========================================================================
 const double TARGET_POINT_X_29 = 19.4; // Coordonnées Point(29) du .geo
 const double TARGET_POINT_Y_29 = 171.0;
 const double TARGET_NODE_TOLERANCE_VISU = 3; // Tolérance pour avertissement dans findClosestNode
 
 // =========================================================================
 //                            Helper Functions
 // =========================================================================
 double funIntegrateArea(double x, double y) { return 1.0; }
 
 // Fonction pour trouver le nœud le plus proche d'un point cible
 int findClosestNode(femNodes *theNodes, double targetX, double targetY) {
     if (!theNodes || theNodes->nNodes == 0) return -1;
     int closestNode = -1;
     double minDistSq = 1e30;
     for (int i = 0; i < theNodes->nNodes; ++i) {
         // Utilise les coordonnées ACTUELLES des noeuds (X, Y)
         double dx = theNodes->X[i] - targetX;
         double dy = theNodes->Y[i] - targetY;
         double distSq = dx * dx + dy * dy;
         if (distSq < minDistSq) {
             minDistSq = distSq;
             closestNode = i;
         }
     }
      // Avertissement si le nœud trouvé est plus loin que la tolérance
     if (sqrt(minDistSq) > TARGET_NODE_TOLERANCE_VISU) {
          fprintf(stderr,"    Warning: Closest node %d (dist %.3e) to target (%.2f, %.2f) is further than tolerance %.3f.\n",
                  closestNode, sqrt(minDistSq), targetX, targetY, TARGET_NODE_TOLERANCE_VISU);
     }
     return closestNode;
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
     printf("    Target point P29(%.2f, %.2f) displacement will be shown in visu.\n", TARGET_POINT_X_29, TARGET_POINT_Y_29);
     printf("\n\n");
 
 
     // =========================================================================
     //                      STAGE 1: Initialization
     // =========================================================================
     printf("STAGE 1: Initializing Geometry and FEM Library...\n");
     geoInitialize();
     femGeo* theGeometry = geoGetGeometry();
     theGeometry->h = 30.0; // Valeur de h pour cette simulation unique
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
     //               STAGE 4: Optional Mesh Saving / Visu Setup
     // =========================================================================
     printf("STAGE 4: Optional Saving and Mesh Visualization Setup...\n");
     printf("  Writing final (fixed) mesh to elasticity_final.txt...\n");
     geoMeshWrite("../data/mesh.txt"); // Vérifiez le chemin
     femNodes *theNodes = theGeometry->theNodes; // Raccourci
     double *meshSizeField = malloc(theNodes->nNodes * sizeof(double));
     if (!meshSizeField) { Error("Failed to allocate meshSizeField"); }
     printf("  Calculating mesh size field for visualization...\n");
     for(int i=0; i < theNodes->nNodes; ++i) { meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]); }
 
     // Bloc de visualisation optionnelle du maillage
     printf("  Mesh visualization block present but disabled (enable by changing 'if(false)').\n");
     if (true) // Activer/déactiver la visualisation du maillage
     {
         printf("  Starting mesh visualization (optional)...\n");
         GLFWwindow* window_mesh = glfemInit("RailUIC60 : Mesh generation (Fixed)");
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
     // Votre bloc femElasticityAddBoundaryCondition original (copiez-le ici)
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
     //femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_Y, 0.0); // Modifié dans votre code
     femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_X, 0.0);
     //femElasticityAddBoundaryCondition(theProblem, "Haut2", DIRICHLET_Y, 0.0); // Modifié dans votre code
     // Charge Neumann (votre code récent)
     double load = -100000;
     femElasticityAddBoundaryCondition(theProblem, "Haut2", NEUMANN_Y, load*0.5);
     femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, load*0.25);
     femElasticityAddBoundaryCondition(theProblem, "Haut1", NEUMANN_Y, load*0.25);
     // Patine1 Diri (votre code récent)
     femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_X, 0.0);
     femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_Y, 0.0);
     
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
     //      ====> STAGE 7bis: Find Target Node & Get Displacement <====
     // =========================================================================
     printf("STAGE 7bis: Finding target node P29 and its displacement...\n");
     // On utilise theNodes qui contient encore les coordonnées ORIGINALES
     int targetNodeIndex29 = findClosestNode(theNodes, TARGET_POINT_X_29, TARGET_POINT_Y_29);
     double ux_target29 = NAN; // Initialise à Not-a-Number
     double uy_target29 = NAN; // Initialise à Not-a-Number
 
     if (targetNodeIndex29 != -1) {
         int dofIndexX = 2 * targetNodeIndex29 + 0;
         int dofIndexY = 2 * targetNodeIndex29 + 1;
         // Vérifier que les indices sont valides pour le vecteur solution
         if (dofIndexX < theProblem->system->size && dofIndexY < theProblem->system->size) {
             ux_target29 = theSoluce[dofIndexX];
             uy_target29 = theSoluce[dofIndexY];
             printf("  Displacement at node %d (closest to P29): Ux = %.6e, Uy = %.6e\n", targetNodeIndex29, ux_target29, uy_target29);
         } else {
             fprintf(stderr, "  Warning: Could not retrieve displacement for target node %d (DOF indices %d, %d out of bounds %d).\n",
                     targetNodeIndex29, dofIndexX, dofIndexY, theProblem->system->size);
             targetNodeIndex29 = -1; // Marquer comme invalide si on ne peut lire la solution
         }
     } else {
         fprintf(stderr, "  Warning: Node closest to P29 (%.2f, %.2f) not found.\n", TARGET_POINT_X_29, TARGET_POINT_Y_29);
     }
 
 
     // =========================================================================
     //         STAGE 8: Post-Processing and Visualization Setup
     // =========================================================================
     printf("STAGE 8: Post-Processing Results for Visualization...\n");
     double deformationFactor = 30000.0; // Vos valeurs
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
     double *originalX = malloc(theNodes->nNodes * sizeof(double));
     double *originalY = malloc(theNodes->nNodes * sizeof(double));
     if (!normDisplacement || !forcesX || !forcesY || !originalX || !originalY) { Error("Failed to allocate visualization fields"); }
 
     printf("  Processing nodal results and applying deformation factor...\n");
     for (int i = 0; i < theNodes->nNodes; i++) {
         originalX[i] = theNodes->X[i]; // Sauvegarde avant modif
         originalY[i] = theNodes->Y[i];
         if (theSoluce && theForces && (2*i + 1) < theProblem->system->size) {
             double ux = theSoluce[2 * i + 0];
             double uy = theSoluce[2 * i + 1];
             theNodes->X[i] += ux * deformationFactor; // Applique déformation pour visu
             theNodes->Y[i] += uy * deformationFactor;
             normDisplacement[i] = sqrt(ux*ux + uy*uy);
             forcesX[i] = theForces[2 * i + 0];
             forcesY[i] = theForces[2 * i + 1];
         } else {
             fprintf(stderr, "  Warning: Error accessing solution/forces for node %d.\n", i);
             normDisplacement[i] = 0.0; forcesX[i] = 0.0; forcesY[i] = 0.0;
             theNodes->X[i] = originalX[i]; // Garder les coordonnées originales si erreur
             theNodes->Y[i] = originalY[i];
         }
     }
 
 
     // =========================================================================
     //                 STAGE 9: Output Numerical Results
     // =========================================================================
     printf("STAGE 9: Outputting Numerical Results...\n");
     // geoMeshWrite reste inchangé...
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
     
     // Afficher aussi le déplacement du point cible dans la console
     if (targetNodeIndex29 != -1) {
        printf("  ==== Target P29 (Node %d) Uy       : %14.7e [m] \n", targetNodeIndex29, uy_target29);
     } else {
        printf("  ==== Target P29 Uy                 : Node not found \n");
     }
 
 
     // =========================================================================
     //                 STAGE 10: Results Visualization
     // =========================================================================
     printf("STAGE 10: Starting Results Visualization (OpenGL)...\n");
     GLFWwindow* window = glfemInit("RailUIC60 : Linear Elasticity Results");
     if (!window) { /* ... gestion erreur ... */ }
     else {
         glfwMakeContextCurrent(window);
         int mode = 1; int domain = 0; int freezingButton = FALSE;
         double t, told = 0.0; char theMessage[256];
 
         do {
            int w, h;
            glfwGetFramebufferSize(window, &w, &h);
            // Attention: glfemReshapeWindows utilise theNodes qui sont maintenant DEFORMES
            glfemReshapeWindows(theNodes, w, h);
            t = glfwGetTime();
 
            // Gestion des clés reste inchangée...
            if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { mode = 0; }
            if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) { mode = 1; }
            if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) { mode = 2; }
            if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) { mode = 3; }
            if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS && !freezingButton) {
               if (theGeometry->nDomains > 0) domain = (domain + 1) % theGeometry->nDomains;
               freezingButton = TRUE; told = t;
            }
            if (t - told > 0.5) { freezingButton = FALSE; }
 
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
 
            // Modification de l'affichage du message
            switch(mode) {
               case 0: // Affichage Domaine
                  if (theGeometry->nDomains > 0 && domain < theGeometry->nDomains && theGeometry->theDomains[domain]) {
                     glfemPlotDomain(theGeometry->theDomains[domain]);
                     // Ajouter info P29 au message du domaine ? Peut-être pas pertinent ici.
                     sprintf(theMessage, "Domain %d: %s (%d edges)", domain, theGeometry->theDomains[domain]->name, theGeometry->theDomains[domain]->nElem);
                  } else { sprintf(theMessage, "No domains or invalid index %d", domain); }
                  break;
               case 1: // Affichage Norme Déplacement
                  glfemPlotField(theGeometry->theElements, normDisplacement);
                  glfemPlotMesh(theGeometry->theElements);
                  // Modifier le message pour inclure Uy de P29
                  if (targetNodeIndex29 != -1) { // Vérifier si le noeud a été trouvé
                     sprintf(theMessage, "Disp Norm (Min:%.2e, Max:%.2e) | P29(Node %d) Uy:%.3e (Def x%.0f)",
                           hMin, hMax, targetNodeIndex29, uy_target29, deformationFactor);
                  } else {
                     sprintf(theMessage, "Disp Norm (Min:%.2e, Max:%.2e) | P29 Node Not Found (Def x%.0f)",
                           hMin, hMax, deformationFactor);
                  }
                  break;
               case 2: // Affichage Force X
                  glfemPlotField(theGeometry->theElements, forcesX);
                  glfemPlotMesh(theGeometry->theElements);
                   // Ajouter info P29 au message des forces X ?
                 if (targetNodeIndex29 != -1) {
                     sprintf(theMessage, "Residual Forces X | P29(Node %d) Uy:%.3e (Def x%.0f)",
                            targetNodeIndex29, uy_target29, deformationFactor);
                 } else {
                     sprintf(theMessage, "Residual Forces X | P29 Node Not Found (Def x%.0f)",
                            deformationFactor);
                 }
                  break;
               case 3: // Affichage Force Y
                  glfemPlotField(theGeometry->theElements, forcesY);
                  glfemPlotMesh(theGeometry->theElements);
                  // Ajouter info P29 au message des forces Y ?
                 if (targetNodeIndex29 != -1) {
                     sprintf(theMessage, "Residual Forces Y | P29(Node %d) Uy:%.3e (Def x%.0f)",
                            targetNodeIndex29, uy_target29, deformationFactor);
                 } else {
                    sprintf(theMessage, "Residual Forces Y | P29 Node Not Found (Def x%.0f)",
                           deformationFactor);
                 }
                  break;
               default:
                   sprintf(theMessage, "Unknown display mode: %d", mode);
                   break;
            }
 
            glColor3f(0.0f, 0.0f, 0.0f);
            glfemMessage(theMessage); // Afficher le message mis à jour
 
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
         // IMPORTANT: Vérifier que theNodes existe toujours et a le bon nombre de noeuds
         if (theNodes && theNodes->nNodes > 0) {
             int expected_nNodes = theNodes->nNodes; // Capturer avant la boucle
              for(int i=0; i < expected_nNodes; ++i) {
                   // Double vérification des indices peut être utile si nNodes a changé ?
                   // Mais normalement, il ne devrait pas changer entre STAGE 8 et STAGE 11.
                   theNodes->X[i] = originalX[i];
                   theNodes->Y[i] = originalY[i];
               }
         } else {
              fprintf(stderr, "  Warning: Could not restore original coordinates (theNodes invalid or empty).\n");
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
     if (theProblem) { // Vérifier s'il n'a pas déjà été libéré (ne devrait pas arriver ici)
         femElasticityFree(theProblem); theProblem = NULL;
     }
     // Libérer la géométrie/maillage et finaliser Gmsh
     geoFinalize();
     // Terminer GLFW
     glfwTerminate();
     printf("Program finished successfully.\n");
     return 0;
 }