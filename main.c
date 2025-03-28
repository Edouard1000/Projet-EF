#include "glfem.h"


// Fonction d'intégration (pour le calcul du poids)
double fun(double x, double y) {
    return 1;
}

int main(void) {
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");

    // Initialisation de la géométrie et de Gmsh
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();

    theGeometry->h = 3*10;
    theGeometry->elementType = FEM_TRIANGLE;
    
    geoMeshGenerate(); // Génère le maillage avec Gmsh
    geoMeshImport();  // Importe le maillage dans les structures de données de fem.c
    
    // Nommer les domaines
    // Domain Bottom
    geoSetDomainName(22, "Bottom0");
    geoSetDomainName(4, "Bottom1");
    geoSetDomainName(5, "Bottom2");
    geoSetDomainName(6, "Bottom3");
    geoSetDomainName(0, "Bottom4");
    geoSetDomainName(1, "Bottom5");
    geoSetDomainName(2, "Bottom6");
    geoSetDomainName(3, "Bottom7");
    geoSetDomainName(17, "Bottom8");

    // Domain Top
    geoSetDomainName(15, "Patine0");
    geoSetDomainName(16, "Patine1");

    // Domain Gauche
    geoSetDomainName(8, "Gauche0");
    geoSetDomainName(10, "Gauche1");
    geoSetDomainName(14, "Gauche2");


    // Domain Droite
    geoSetDomainName(7, "Droite0");
    geoSetDomainName(9, "Droite1");
    geoSetDomainName(11, "Droite2");

    // Domain Haut
    geoSetDomainName(12, "Haut0");
    geoSetDomainName(21, "Haut1");
    geoSetDomainName(20, "Haut2");
    geoSetDomainName(19, "Haut3");
    geoSetDomainName(18, "Haut4");
    geoSetDomainName(13, "Haut5");




    geoMeshWrite("../data/elasticity.txt"); // Optionnel: Sauvegarde du maillage
 //
 //  -3- Champ de la taille de r�f�rence du maillage
 //

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes*sizeof(double));
    femNodes *theNodes = theGeometry->theNodes;
    for(int i=0; i < theNodes->nNodes; ++i)
        meshSizeField[i] = geoSize(theNodes->X[i], theNodes->Y[i]);
 //
 //  -4- Visualisation du maillage
 //  
     
    int mode = 1; // Change mode by pressing "j", "k", "l"
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[256];
    double pos[2] = {20,460};


    GLFWwindow* window = glfemInit("EPL1110 : Mesh generation ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;


        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
    //    glfemChangeState(&mode, theMeshes->nMesh);
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}

        
        if (t-told > 0.5) {freezingButton = FALSE; }
            
        
        
        
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);

            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);

            
            
            }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            
            
            
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);

            
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
            }
            
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
            glfwWindowShouldClose(window) != 1 );

    // Paramètres du problème d'élasticité
    double E   = 211e9;     // Module de Young
    double nu  = 0.3;      // Coefficient de Poisson
    double rho = 7850;      // Masse volumique
    double g   = 9.81;     // Accélération gravitationnelle


    // Création du problème
    femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);
    printf("femProblem size : %d \n",theProblem->system->size);
    // Ajout des conditions aux limites (Dirichlet et Neumann)

    // Add boundary conditions for Bottom domains
    femElasticityAddBoundaryCondition(theProblem, "Bottom0", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom2", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom3", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom4", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom5", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom6", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom7", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Bottom8", DIRICHLET_Y, 0.0);

    // Add boundary conditions for Gauche domains
    femElasticityAddBoundaryCondition(theProblem, "Gauche0", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Gauche1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Gauche2", DIRICHLET_X, 0.0);

    // Add boundary conditions for Droite domains
    femElasticityAddBoundaryCondition(theProblem, "Droite0", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite1", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Droite2", DIRICHLET_X, 0.0);

    // Add boundary conditions for Haut domains
    femElasticityAddBoundaryCondition(theProblem, "Haut0", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut1", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut2", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut3", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut4", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Haut5", DIRICHLET_Y, 0.0);



    femElasticityAddBoundaryCondition(theProblem, "Patine0", NEUMANN_Y, -1e4);
    femElasticityAddBoundaryCondition(theProblem, "Patine1", NEUMANN_Y, -1e4);
    


    femElasticityPrint(theProblem); // Affiche les paramètres du problème


    // Résolution du problème
    double *theSoluce = femElasticitySolve(theProblem);

    // Calcul des forces résiduelles
    double *theForces = femElasticityForces(theProblem);

     // Intégration de la fonction 'fun' (pour le calcul du poids total)
    double area = femElasticityIntegrate(theProblem, fun);

    // Facteur d'amplification pour la visualisation des déplacements
    //femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e4; // A ajuster selon votre géométrie/déplacements

    // Allocation de mémoire pour les champs à afficher
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    
    
    for (int i = 0; i < theNodes->nNodes; i++) {
        // Applique la déformation pour la visualisation
        theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
        theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
        // Calcule la norme du déplacement
        normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] +
                                   theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
        // Stocke les forces résiduelles
        forcesX[i] = theForces[2 * i + 0];
        forcesY[i] = theForces[2 * i + 1];
    }

    
    // Affichage des déplacements min et max
    double hMin = femMin(normDisplacement, theNodes->nNodes);
    double hMax = femMax(normDisplacement, theNodes->nNodes);
    printf(" ==== Minimum displacement          : %14.7e [m] \n", hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n", hMax);

    // Calcul et affichage des forces globales
    double theGlobalForce[2] = {0, 0};
    for (int i = 0; i < theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2 * i + 0];
        theGlobalForce[1] += theForces[2 * i + 1];
    }
    printf(" ==== Global horizontal force       : %14.7e [N] \n", theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n", theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

     // Initialisation de la fenêtre GLFW
    // GLFWwindow* window = glfemInit("EPL1110 : Linear Elasticity");
    // glfwMakeContextCurrent(window);

    // int mode = 1; // Mode d'affichage
    // int domain = 0; // Domaine courant
    // int freezingButton = FALSE;
    // double t, told = 0;
    // char theMessage[256];

    // // Boucle principale de l'application (boucle d'événements)
    // do {
    //     int w, h;
    //     glfwGetFramebufferSize(window, &w, &h); // Récupère la taille de la fenêtre
    //     glfemReshapeWindows(theGeometry->theNodes, w, h); // Adapte la vue

    //     t = glfwGetTime(); // Temps courant pour les animations

    //     // Gestion des touches du clavier pour changer le mode d'affichage
    //     if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
    //     if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
    //     if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 2; }
    //     if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 3; }
    //     if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) {
    //         domain++;
    //         freezingButton = TRUE;
    //         told = t;
    //     }
    //     if (t - told > 0.5) {
    //         freezingButton = FALSE;
    //     }

    //      // Affichage en fonction du mode
    //     if (mode == 0) {
    //         domain = domain % theGeometry->nDomains;  // Cycle sur les domaines
    //         glfemPlotDomain(theGeometry->theDomains[domain]);  // Affiche le domaine courant
    //         sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
    //         glColor3f(1.0, 0.0, 0.0);
    //         glfemMessage(theMessage); // Affiche le nom du domaine
    //     }
    //     if (mode == 1) {
    //          // Affiche le champ de déplacement et le maillage déformé
    //         glfemPlotField(theGeometry->theElements, normDisplacement);
    //         glfemPlotMesh(theGeometry->theElements);
    //         sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
    //         glColor3f(1.0, 0.0, 0.0);
    //         glfemMessage(theMessage);  // Affiche le nombre d'éléments
    //     }
    //     if (mode == 2) {
    //          // Affiche les forces résiduelles horizontales et le maillage
    //         glfemPlotField(theGeometry->theElements, forcesX);
    //         glfemPlotMesh(theGeometry->theElements);
    //         sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
    //         glColor3f(1.0, 0.0, 0.0);
    //         glfemMessage(theMessage); // Affiche le nombre d'éléments
    //     }
    //     if (mode == 3) {
    //         // Affiche les forces résiduelles verticales et le maillage
    //         glfemPlotField(theGeometry->theElements, forcesY);
    //         glfemPlotMesh(theGeometry->theElements);
    //         sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
    //         glColor3f(1.0, 0.0, 0.0);
    //         glfemMessage(theMessage);  // Affiche le nombre d'éléments
    //     }
    //     glfwSwapBuffers(window); // Swap des buffers (double buffering)
    //     glfwPollEvents();       // Traite les événements (clavier, souris, etc.)

    // } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
    //          glfwWindowShouldClose(window) != 1);

    // Libération de la mémoire (très important !)
    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem); // Libère *toutes* les structures du problème
    geoFinalize();      // Finalise Gmsh et libère la géométrie
    glfwTerminate();    // Termine GLFW

    exit(EXIT_SUCCESS); // Fin du programme
    return 0;
}