/**
 * @file main.c
 * @brief Programme principal pour la simulation par éléments finis d'un rail UIC60
 * @authors Edouard Meurant & Théodore Moulaert
 *
 * Ce programme effectue les étapes suivantes :
 * 1. Initialisation de la bibliothèque FEM et Gmsh.
 * 2. Génération, importation et correction du maillage (UIC60.geo par défaut).
 * 3. Attribution de noms aux frontières (domaines).
 * 4. Visualisation optionnelle du maillage (commentée par défaut).
 * 5. Configuration du problème d'élasticité (Planar Strain OU Axisymétrique).
 * 6. Application des conditions aux limites (Dirichlet et Neumann) adaptées au cas choisi.
 * 7. Résolution du système linéaire par éléments finis.
 * 8. Post-traitement des résultats (calculs, préparation visualisation).
 * 9. Sauvegarde des déplacements nodaux dans data/elasticity.txt.
 * 10. Visualisation interactive des résultats via OpenGL.
 * 11. Nettoyage et libération de la mémoire.
 */

 #include "glfem.h"    // Pour la visualisation OpenGL
 #include <stdlib.h>   // Pour malloc, free, exit
 #include <stdio.h>    // Pour printf, fprintf, fopen, fclose
 #include <math.h>     // Pour sqrt, M_PI (si _USE_MATH_DEFINES est activé dans homework.c)
 #include <stdbool.h>  // Pour le type bool et true/false
 


 // =========================================================================
 //                            Target Point Definition
 // =========================================================================
 const double TARGET_POINT_X_29 = 19.4;         // Coordonnée x Point(29) du .geo
 const double TARGET_POINT_Y_29 = 171.0;        // Coordonnée y Point(29) du .geo
 const double TARGET_NODE_TOLERANCE_VISU = 2;   // Tolérance pour avertissement dans findClosestNode 
 

 // =========================================================================
 //                            Helper Function
 // =========================================================================
 
 /**
  * @brief Fonction constante pour l'intégration de l'aire (ou volume en axisym).
  * Attention: le résultat n'est physiquement correct (aire/volume) que si
  * la fonction d'intégration dans fem.c/homework.c est adaptée.
  */
 double funIntegrateArea(double x, double y) {
     return 1.0;
 }



 /**
  * @brief Trouve le nœud le plus proche d'un point cible (targetX, targetY).
  * @param theNodes Pointeur vers la structure de nœuds.
  * @param targetX Coordonnée x du point cible.
  * @param targetY Coordonnée y du point cible.
  * @return L'indice du nœud le plus proche ou -1 si aucun nœud trouvé.
  */
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
 
     // --- User Instructions & Welcome Message ---
     printf("\n\n");
     printf("    FEM Simulation - Rail UIC60 / Axisymmetric Test\n");
     printf("    ----------------------------------------------\n");
     printf("    OpenGL Visualization Keys:\n");
     printf("      V : Displacement norm field (on deformed mesh)\n");
     printf("      D : Highlight domains (boundaries)\n");
     printf("      X : Horizontal (or Radial) residual forces field\n");
     printf("      Y : Vertical (or Axial) residual forces field\n");
     printf("      N : Cycle through highlighted domains\n");
     printf("      ESC : Exit visualization window\n");
     printf("\n\n");
 
 
     // =========================================================================
     //                      STAGE 1: Initialization
     // =========================================================================
     printf("STAGE 1: Initializing Geometry and FEM Library...\n");
 
     geoInitialize(); // Initialise Gmsh et la structure theGeometry
     femGeo* theGeometry = geoGetGeometry();
 
     // Paramètres de maillage
     theGeometry->h = 30.0;                     // Taille de maille de référence
     theGeometry->elementType = FEM_TRIANGLE;   // Type d'élément par défaut
 
     // =========================================================================
     //                 STAGE 2: Mesh Generation and Processing
     // =========================================================================
     printf("STAGE 2: Generating, Importing, and Fixing Mesh...\n");
     printf("  Using geometry file: ../Projet-EF/UIC60.geo (default)\n");
     printf("  NOTE: For AXISYMMETRIC case, this geometry might not be physically meaningful.\n");
     // Génération Gmsh (utilise UIC60.geo et la callback geoSize de homework.c)
     geoMeshGenerate();
     // Importation du maillage depuis Gmsh vers la structure femGeo
     geoMeshImport();
     // Correction du maillage (suppression des noeuds non utilisés par les éléments 2D)
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

     // Ces noms sont utilisés pour appliquer les Conditions aux Limites
     if (theGeometry->nDomains > 0)  geoSetDomainName(0, "Bottom4");   else printf("  Warning: Skipping domain index 0 (nDomains=%d)\n", theGeometry->nDomains);
     if (theGeometry->nDomains > 1)  geoSetDomainName(1, "Bottom5");
     if (theGeometry->nDomains > 2)  geoSetDomainName(2, "Bottom6");
     if (theGeometry->nDomains > 3)  geoSetDomainName(3, "Bottom7");
     if (theGeometry->nDomains > 4)  geoSetDomainName(4, "Bottom1");
     if (theGeometry->nDomains > 5)  geoSetDomainName(5, "Bottom2");
     if (theGeometry->nDomains > 6)  geoSetDomainName(6, "Bottom3");
     if (theGeometry->nDomains > 7)  geoSetDomainName(7, "Droite0");    // Potentiellement surface externe r=max en AXISYM
     if (theGeometry->nDomains > 8)  geoSetDomainName(8, "Gauche0");    // Potentiellement axe r=0 en AXISYM
     if (theGeometry->nDomains > 9)  geoSetDomainName(9, "Droite1");
     if (theGeometry->nDomains > 10) geoSetDomainName(10, "Gauche1");   // Potentiellement axe r=0 en AXISYM
     if (theGeometry->nDomains > 11) geoSetDomainName(11, "Droite2");
     if (theGeometry->nDomains > 12) geoSetDomainName(12, "Haut0");     // Surface supérieure (charge en planaire)
     if (theGeometry->nDomains > 13) geoSetDomainName(13, "Haut5");
     if (theGeometry->nDomains > 14) geoSetDomainName(14, "Gauche2");   // Potentiellement axe r=0 en AXISYM
     if (theGeometry->nDomains > 15) geoSetDomainName(15, "Patine0");   // Surface supérieure (charge en planaire)
     if (theGeometry->nDomains > 16) geoSetDomainName(16, "Patine1");
     if (theGeometry->nDomains > 17) geoSetDomainName(17, "Bottom8");
     if (theGeometry->nDomains > 18) geoSetDomainName(18, "Haut4");
     if (theGeometry->nDomains > 19) geoSetDomainName(19, "Haut3");
     if (theGeometry->nDomains > 20) geoSetDomainName(20, "Haut2");     // Surface supérieure (charge en planaire)
     if (theGeometry->nDomains > 21) geoSetDomainName(21, "Haut1");     // Surface supérieure (charge en planaire)
     if (theGeometry->nDomains > 22) geoSetDomainName(22, "Bottom0"); else if (theGeometry->nDomains <= 22) printf("  Warning: Skipping domain index 22\n");
     printf("  Domain naming attempt completed.\n");
 
 
     // =========================================================================
     //               STAGE 4: Optional Mesh Saving / Visu Setup
     // =========================================================================
     printf("STAGE 4: Optional Saving and Mesh Visualization Setup...\n");
 
     // Sauvegarde du maillage final (après correction) pour inspection externe
     printf("  Writing final (fixed) mesh to data/mesh.txt...\n");
     geoMeshWrite("../data/mesh.txt");
 
     // Calcul du champ de taille de référence pour visualisation
     femNodes *theNodes = theGeometry->theNodes;
     double *meshSizeField = malloc(theNodes->nNodes * sizeof(double));
     if (!meshSizeField) { Error("Failed to allocate meshSizeField"); }
     printf("  Calculating mesh size field for visualization...\n");
     for(int i=0; i < theNodes->nNodes; ++i) {
         meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
     }
 
     // ---Bloc de visualisation optionnelle du maillage---
     printf("  Mesh visualization block present but disabled (enable by changing 'if(false)').\n");
     if (false) // Activer/déactiver la visualisation du maillage
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
 
     double E   = 211e9; // Module d'Young [Pa]
     double nu  = 0.3;   // Coefficient de Poisson [-]
     double rho = 10000; // Masse volumique [kg/m3]
     double g   = 9.81;  // Accélération due à la gravité [m/s2] (mettre 0 si non désirée)
 
     // ---- CHOIX DU TYPE DE PROBLEME ----
     femElasticCase problemCase = PLANAR_STRAIN;    // Cas par défaut: Rail UIC60
     //femElasticCase problemCase = AXISYM;         // Décommenter pour tester AXISYM
 
     printf("  Creating femProblem for case: %s\n",
            (problemCase == PLANAR_STRAIN) ? "PLANAR_STRAIN" :
            (problemCase == PLANAR_STRESS) ? "PLANAR_STRESS" : "AXISYM");
 
     femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, problemCase);
     if (!theProblem) {
         fprintf(stderr, "Error: femElasticityCreate failed!\n");
         free(meshSizeField); // Libérer mémoire allouée avant exit
         geoFinalize();
         exit(EXIT_FAILURE);
     }
     printf("  femProblem created. System size (DOFs): %d \n", theProblem->system->size);
 
 
     // =========================================================================
     //                 STAGE 6: Apply Boundary Conditions
     // =========================================================================
     printf("STAGE 6: Applying Boundary Conditions...\n");
 
     if (problemCase == AXISYM) {
         printf("  Applying AXISYMMETRIC Boundary Conditions...\n");
         // Interprétation: X->Radial (r), Y->Axial (z)
         // 'Gauche*' = Axe de symétrie (r=0)
         // 'Bottom*' = Base fixe (z=0)
 
         // --- Conditions aux limites Dirichlet AXISYMETRIQUES ---
 
         // Axe de Symétrie (r=0) : Déplacement radial Ur = 0 (DIRICHLET_X)
         // Correspond aux domaines 'Gauche*' de la géométrie UIC60
         // Commentaire: Le déplacement axial Uz (DIRICHLET_Y) est laissé libre sur l'axe.
         femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_X, 0.0);
         femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_X, 0.0);
         femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);
 
         // Base Fixe (z=0) : Déplacement axial Uz = 0 (DIRICHLET_Y)
         //                     ET Déplacement radial Ur = 0 (DIRICHLET_X)
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
        femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_Y, 0.0);
        
        // Droite
        femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_Y, 0.0);
        
        // Haut
        femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_X, 0.0);
        
        // Patine
        femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_Y, 0.0);
        
        // ---Charge Neumann---
        double load = -6250000; // Valeur de la charge appliquée [N/m]

        femElasticityAddBoundaryCondition(theProblem, "Haut2", NEUMANN_Y, load*0.5);
        femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, load*0.25);
        femElasticityAddBoundaryCondition(theProblem, "Haut1", NEUMANN_Y, load*0.25);

         
     } else { // Cas PLANAR_STRAIN
        printf("  Applying PLANAR Boundary Conditions (for UIC60 Rail)...\n");

        // ---Conditions aux limites Dirichlet---
     
        // Bas
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
        femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_Y, 0.0);
        
        // Droite
        femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_Y, 0.0);
        
        // Haut
        femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_X, 0.0);
        
        // Patine
        femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_X, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Patine1", DIRICHLET_Y, 0.0);
        
        // ---Charge Neumann---
        double load = -6250000; // Valeur de la charge appliquée [N/m]

        femElasticityAddBoundaryCondition(theProblem, "Haut2", NEUMANN_Y, load*0.5);
        femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, load*0.25);
        femElasticityAddBoundaryCondition(theProblem, "Haut1", NEUMANN_Y, load*0.25);
     }
 
     // Affichage du résumé du problème après application des CL
     femElasticityPrint(theProblem);
 
 
     // =========================================================================
     //                 STAGE 7: Solve the FEM Problem
     // =========================================================================
     printf("STAGE 7: Assembling and Solving the Linear System...\n");
     
     femSolver *gaussSolver = malloc(sizeof(femSolver));
     gaussSolver->type = FEM_FULL;
     gaussSolver->solver = theProblem->system;

     double *theSoluce = femElasticitySolve(theProblem);
     if (!theSoluce) { // Vérifier si la solution a échoué (e.g. pivot nul dans l'élimination)
         fprintf(stderr, "Error: FEM solver failed (femElasticitySolve returned NULL). Check BCs (rigid body motion?) or matrix assembly.\n");
         femElasticityFree(theProblem);
         free(meshSizeField);
         geoFinalize();
         exit(EXIT_FAILURE);
     }
     printf("  System solved.\n");
 
     printf("  Calculating residual forces...\n");
     double *theForces = femElasticityForces(theProblem);
     printf("  Residual forces calculated.\n");
 
     printf("  Calculating area/volume integral (f=1)...\n");
     double area = femElasticityIntegrate(theProblem, funIntegrateArea); // Intégrale de 1 sur le domaine 2D
     if (problemCase == AXISYM) {
         printf("  NOTE (Axisymmetric): The value %.4e represents Integral(2*pi*r dr dz).\n", area);
         printf("                     This corresponds to the volume of the generated 3D shape.\n");
     } else {
         printf("  NOTE (Planar): The value %.4e represents the area of the 2D section [m^2].\n", area);
     }

     


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
 
     // Facteur d'amplification pour la visualisation des déformations
     double deformationFactor = 1.0; // Valeur par défaut
     if (problemCase == PLANAR_STRAIN || problemCase == PLANAR_STRESS) {
         deformationFactor = 10000.0; // Amplification importante pour voir déformation du rail
         printf("  Using deformation factor for visualization: %.1f (Planar Case)\n", deformationFactor);
     } else if (problemCase == AXISYM) {
         deformationFactor = 1000.0;
          printf("  Using deformation factor for visualization: %.1f (Axisymmetric Case)\n", deformationFactor);
     }
 
     // Allocation mémoire pour les champs nodaux de visualisation
     double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
     double *forcesX = malloc(theNodes->nNodes * sizeof(double));
     double *forcesY = malloc(theNodes->nNodes * sizeof(double));
     double *originalX = malloc(theNodes->nNodes * sizeof(double)); // Sauvegarde coords originales
     double *originalY = malloc(theNodes->nNodes * sizeof(double));
      if (!normDisplacement || !forcesX || !forcesY || !originalX || !originalY) {
          // Libérer ce qui a pu être alloué
          free(normDisplacement); free(forcesX); free(forcesY); free(originalX);
          Error("Failed to allocate visualization/backup fields");
      }
 
     // Calcul des champs et application du facteur de déformation pour la visu
     printf("  Processing nodal results and applying deformation factor...\n");
     for (int i = 0; i < theNodes->nNodes; i++) {
         originalX[i] = theNodes->X[i]; // Sauvegarde avant modification
         originalY[i] = theNodes->Y[i];
         if (theSoluce && theForces && (2*i + 1) < theProblem->system->size) {
             double ux = theSoluce[2 * i + 0]; // Déplacement X ou R
             double uy = theSoluce[2 * i + 1]; // Déplacement Y ou Z
             theNodes->X[i] += ux * deformationFactor; // Modifie pour visu
             theNodes->Y[i] += uy * deformationFactor; // Modifie pour visu
             normDisplacement[i] = sqrt(ux*ux + uy*uy);
             forcesX[i] = theForces[2 * i + 0]; // Force X ou R
             forcesY[i] = theForces[2 * i + 1]; // Force Y ou Z
         } else {
              fprintf(stderr, "  Warning: Error accessing solution/forces for node %d.\n", i);
              normDisplacement[i] = 0.0; forcesX[i] = 0.0; forcesY[i] = 0.0;
              // Ne pas appliquer de déformation si données invalides
              theNodes->X[i] = originalX[i];
              theNodes->Y[i] = originalY[i];
         }
     }
 
 
     // =========================================================================
     //                 STAGE 9: Output Numerical Results
     // =========================================================================
     printf("STAGE 9: Outputting Numerical Results...\n");
 
     // Sauvegarde des déplacements nodaux dans un fichier texte
     FILE *solFile = fopen("../data/elasticity.txt", "w");
     if (solFile) {
         fprintf(solFile, "FEM Elasticity Results\n");
         fprintf(solFile, "Problem Type: %s\n",
                 (problemCase == PLANAR_STRAIN) ? "PLANAR_STRAIN" :
                 (problemCase == PLANAR_STRESS) ? "PLANAR_STRESS" : "AXISYM");
         fprintf(solFile, "Number of nodes: %d\n", theGeometry->theNodes->nNodes);
         fprintf(solFile, "Number of DOFs: %d\n", theProblem->system->size);
         fprintf(solFile, "Node#\tOrigX(R)\tOrigY(Z)\tDispX(R)\tDispY(Z)\n");
         for (int i=0; i < theGeometry->theNodes->nNodes; ++i) {
             fprintf(solFile, "%d\t%.6e\t%.6e\t%.6e\t%.6e\n",
                     i,
                     originalX[i], // Coordonnées originales
                     originalY[i],
                     (2*i < theProblem->system->size) ? theSoluce[2*i] : 0.0,
                     (2*i+1 < theProblem->system->size) ? theSoluce[2*i+1] : 0.0);
         }
         fclose(solFile);
         printf("  Solution (displacements) written to data/elasticity.txt\n");
     } else {
         printf("  Warning: Could not open data/elasticity.txt for writing solution.\n");
     }
 
     // Calcul et affichage des déplacements min/max
     double dMin = femMin(normDisplacement, theNodes->nNodes);
     double dMax = femMax(normDisplacement, theNodes->nNodes);
     printf("  ==== Minimum displacement norm     : %14.7e [m] \n", dMin);
     printf("  ==== Maximum displacement norm     : %14.7e [m] \n", dMax);
 
     // Calcul et affichage des forces résiduelles globales
     double theGlobalForce[2] = {0.0, 0.0};
     if (theForces) {
         for (int i = 0; i < theNodes->nNodes; i++) {
             // Vérifier si l'indice est valide (pourrait ne pas l'être si le système est plus petit que 2*nNodes)
              if ((2*i + 1) < theProblem->system->size) {
                 theGlobalForce[0] += theForces[2 * i + 0]; // Somme forces X ou R
                 theGlobalForce[1] += theForces[2 * i + 1]; // Somme forces Y ou Z
             }
         }
     }
     printf("  ==== Sum of Global Residuals X (or R): %14.7e [N] \n", theGlobalForce[0]);
     printf("  ==== Sum of Global Residuals Y (or Z): %14.7e [N] \n", theGlobalForce[1]);
     // Affichage du poids estimé (sens physique seulement en planaire)
     if (problemCase != AXISYM) {
        printf("  ==== Estimated Weight (rho*g*Area*1m): %14.7e [N] \n", area * rho * g);
     } else {
         printf("  ==== Estimated Weight (rho*g*Volume): %14.7e [N] \n", area * rho * g);
     }

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
 
     GLFWwindow* window = glfemInit("EPL1110 : Linear Elasticity Results");
     if (!window) {
         fprintf(stderr, "  Error: Failed to initialize GLFW window for results visualization.\n");
         // Pas besoin de quitter, le programme peut continuer sans visu
     } else {
         glfwMakeContextCurrent(window);
         int mode = 1; // Mode d'affichage (0=Domains, 1=DispNorm, 2=ForcesX/R, 3=ForcesY/Z)
         int domain = 0; // Index du domaine affiché en mode 0
         int freezingButton = FALSE;
         double t, told = 0.0;
         char theMessage[256];
 
         do {
             int w, h;
             glfwGetFramebufferSize(window, &w, &h);
             // Utilise les coordonnées nodales modifiées par deformationFactor pour le reshape
             glfemReshapeWindows(theGeometry->theNodes, w, h);
             t = glfwGetTime();
 
             // Gestion des touches (identique à l'original)
             if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) { mode = 0; }
             if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) { mode = 1; }
             if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) { mode = 2; }
             if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) { mode = 3; }
             if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) { mode = 4; }
             if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS && freezingButton == FALSE) {
                 if (theGeometry->nDomains > 0) domain = (domain + 1) % theGeometry->nDomains;
                 freezingButton = TRUE; told = t;
             }
             if (t - told > 0.5) { freezingButton = FALSE; }
 
             glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
             glClear(GL_COLOR_BUFFER_BIT);
 
             // Dessin en fonction du mode
             switch(mode) {
                 case 0: // Affichage des domaines
                     if (theGeometry->nDomains > 0 && domain < theGeometry->nDomains && theGeometry->theDomains[domain]) {
                         glfemPlotDomain(theGeometry->theDomains[domain]); // Utilise les coordonnées déformées
                         sprintf(theMessage, "Domain %d: %s (%d edges)", domain, theGeometry->theDomains[domain]->name, theGeometry->theDomains[domain]->nElem);
                     } else { sprintf(theMessage, "No domains or invalid index %d", domain); }
                     break;
                 case 1: // Norme du déplacement
                     glfemPlotField(theGeometry->theElements, normDisplacement);
                     glfemPlotMesh(theGeometry->theElements); // Maillage déformé
                     sprintf(theMessage, "Displacement Norm (Min: %.2e, Max: %.2e) (Deformed Plot x%.0f)", dMin, dMax, deformationFactor);
                     break;
                 case 2: // Forces X ou R
                     glfemPlotField(theGeometry->theElements, forcesX);
                     glfemPlotMesh(theGeometry->theElements); // Maillage déformé
                     sprintf(theMessage, "Residual Forces X (or R) (Deformed Plot x%.0f)", deformationFactor);
                     break;
                 case 3: // Forces Y ou Z
                     glfemPlotField(theGeometry->theElements, forcesY);
                     glfemPlotMesh(theGeometry->theElements); // Maillage déformé
                     sprintf(theMessage, "Residual Forces Y (or Z) (Deformed Plot x%.0f)", deformationFactor);
                     break;
                 case 4: // Affichage spy-matrice
                     glColor3f(1.0,0.0,0.0);
                     glfemPlotSolver(gaussSolver,theGeometry->theNodes->nNodes,w,h);
                   break;
                 default:
                      sprintf(theMessage, "Unknown display mode: %d", mode);
                      break;
             }
 
             glColor3f(0.0f, 0.0f, 0.0f); // Couleur du texte
             glfemMessage(theMessage); // Affichage du message
 
             glfwSwapBuffers(window); // Afficher le rendu
             glfwPollEvents(); // Traiter les événements (clavier, souris, fermeture)
 
         } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && // Quitter avec ESC
                  glfwWindowShouldClose(window) != 1); // Quitter si fenêtre fermée
 
         glfwDestroyWindow(window); // Fermer la fenêtre OpenGL
         printf("  Results visualization finished.\n");
     }
 
 
     // =========================================================================
     //                 STAGE 11: Cleanup and Finalization
     // =========================================================================
     printf("STAGE 11: Cleaning Up Memory and Libraries...\n");
 
     // Restaurer les coordonnées nodales originales (IMPORTANT!)
     printf("  Restoring original node coordinates...\n");
     if (originalX && originalY) {
         for(int i=0; i < theNodes->nNodes; ++i) {
             theNodes->X[i] = originalX[i];
             theNodes->Y[i] = originalY[i];
         }
         free(originalX); originalX = NULL;
         free(originalY); originalY = NULL;
     } else {
          printf("  Warning: Could not restore original coordinates (backup arrays missing).\n");
     }
 
     // Libérer la mémoire allouée dans main()
     printf("  Freeing allocated memory in main...\n");
     free(meshSizeField);    meshSizeField = NULL;
     free(normDisplacement); normDisplacement = NULL;
     free(forcesX);          forcesX = NULL;
     free(forcesY);          forcesY = NULL;
     // originalX/Y déjà libérés s'ils existaient
 
     // Libérer les structures du problème EF (inclut système, solution, résidus, etc.)
     femElasticityFree(theProblem); theProblem = NULL;
 
     // Libérer la géométrie/maillage et finaliser Gmsh
     geoFinalize();
 
     // Terminer GLFW si la visu a été initialisée
     if (window) { // Vérifier si window a été initialisé (même si échec création)
       glfwTerminate();
     }
 
 
     printf("\nProgram finished successfully.\n");
     return 0; // Succès
 }