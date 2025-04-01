/**
 * @file homework.c
 * @brief Fonctions spécifiques au projet étudiant pour la simulation FEM.
 * @authors Edouard Meurant & Théodore Moulaert

 *
 * Ce fichier contient les implémentations requises pour le projet, notamment :
 * - La fonction de callback pour la taille du maillage (`geoSize`).
 * - Les fonctions d'assemblage de la matrice de rigidité et du vecteur force (`femElasticityAssembleElements`, `femElasticityAssembleNeumann`).
 * - La fonction de résolution du système linéaire (`femElasticitySolve`).
 * - La fonction de calcul des forces résiduelles (`femElasticityForces`).
 * - La fonction d'appel à Gmsh pour la génération du maillage (`geoMeshGenerate`).
 * - Des fonctions utilitaires (interpolation, affichage matrice, vérification symétrie).
 *
 * Il dépend de "fem.h" pour les structures de données et les fonctions de base.
 */

 #include "fem.h"
 #include <math.h>
 #include <stdbool.h> 



 
 // =========================================================================
 //                 Section 1: Geometry and Meshing Functions
 // =========================================================================
 
 /**
  * @brief Interpolation cubique d'Hermite entre h0 et h_star sur une distance d_star.
  *
  * Utilisée par `geoSize` pour assurer une transition douce de la taille de maille.
  *
  * @param d Distance actuelle du point au centre de raffinement.
  * @param d_star Distance d'influence du raffinement.
  * @param h0 Taille de maille cible au centre (pour d=0).
  * @param h_star Taille de maille à l'extérieur de la zone d'influence (pour d>=d_star).
  * @return double Taille de maille interpolée.
  */
 double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
     if (d >= d_star) return h_star;
     double t = d / d_star;
     // Polynôme d'Hermite H1(t) * h0 + H0(t) * h_star, avec H1(t)=3t^2-2t^3, H0(t)=1-H1(t)
     return h0 + (h_star - h0) * (3.0 * t * t - 2.0 * t * t * t);
 }
 
 /**
  * @brief Fonction de callback pour Gmsh définissant la taille de maille souhaitée en (x, y).
  *
  * Cette fonction est appelée par Gmsh lors de la génération du maillage.
  * Elle permet de spécifier une taille de maille variable dans le domaine,
  * typiquement plus fine dans les zones d'intérêt (ici, le haut du rail).
  *
  * @param x Coordonnée X du point où la taille est demandée.
  * @param y Coordonnée Y du point où la taille est demandée.
  * @return double Taille de maille souhaitée au point (x,y).
  */

double pointToSegmentDistance(double x, double y, double x1, double y1, double x2, double y2) {
    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = (len_sq != 0) ? dot / len_sq : -1;

    double xx, yy;

    if (param < 0) {
        xx = x1;
        yy = y1;
    } else if (param > 1) {
        xx = x2;
        yy = y2;
    } else {
        xx = x1 + param * C;
        yy = y1 + param * D;
    }

    double dx = x - xx;
    double dy = y - yy;
    return sqrt(dx * dx + dy * dy);
}

double geoSize(double x, double y) {
    femGeo* theGeometry = geoGetGeometry();
    double h_base = theGeometry->h; 

    // Coordonnées des points d'intérêt pour le raffinement sur le haut du rail UIC60
    const double yTop = 172.0;
    const double X1   = 37.85;
    const double X2   = -X1;
    const double X3   = 33.8;
    const double Y3   = 164.3;
    const double X4   = -X3;
    
    // Segments d'intérêt
    const double xSegStart = 8.25,  ySegStart = 31.5,  xSegEnd = 8.25,  ySegEnd = 121.0;
    const double xSeg2Start = -xSegStart, ySeg2Start = ySegStart, xSeg2End = -xSegEnd, ySeg2End = ySegEnd;

    // Calcul des distances aux zones d'intérêt
    double d_top  = fabs(yTop - y);                                                        // Distance verticale au sommet
    double d_pt1  = fabs(X1 - x);                                                          // Distance horizontale au point 1
    double d_pt2  = fabs(X2 - x);                                                          // Distance horizontale au point 2
    double d_pt3  = sqrt(pow(X3 - x, 2) + pow(Y3 - y, 2));                                 // Distance euclidienne au point 3
    double d_pt4  = sqrt(pow(X4 - x, 2) + pow(Y3 - y, 2));                                 // Distance euclidienne au point 4
    double dSeg   = pointToSegmentDistance(x, y, xSegStart, ySegStart, xSegEnd, ySegEnd);  // Distance au segment 1
    double dSeg2  = pointToSegmentDistance(x, y, xSeg2Start, ySeg2Start, xSeg2End, ySeg2End); // Distance au segment 2

    // Paramètres du raffinement local
    double h_refined = h_base * 0.1;       // Taille de maille cible très fine dans la zone d'intérêt
    double influence_radius = 20.0;        // Rayon d'influence autour des points pour le raffinement

    // Calcul du h dans la zone supérieure
    double h_upper = h_base;
    if (y > 140.0) {
        double min_dist = fmin(fmin(fmin(fmin(d_top, d_pt1), d_pt2), d_pt3), d_pt4);
        h_upper = hermiteInterpolation(min_dist, influence_radius, h_refined, h_base);
    }

    // Raffinement général en fonction des deux segments
    double h_segment  = hermiteInterpolation(dSeg,  25, h_refined, h_base);
    double h_segment2 = hermiteInterpolation(dSeg2, 25, h_refined, h_base);

    // Retour du minimum des influences (haut, segment 1 et segment 2)
    return fmin(h_upper, fmin(h_segment, h_segment2));
}




 
 /**
  * @brief Lance la génération du maillage 2D via l'API Gmsh.
  *
  * Ouvre le fichier .geo spécifié, configure Gmsh avec le callback `geoSize`,
  * et demande à Gmsh de générer le maillage surfacique. Le maillage résultant
  * reste en mémoire dans Gmsh, prêt à être importé via `geoMeshImport()`.
  */
 void geoMeshGenerate() {
     int ierr;
     femGeo* theGeometry = geoGetGeometry(); // Pas utilisé directement ici, mais assure que la structure existe
 
     const char* geoFilePath = "../Projet-EF/UIC60.geo";
 
     printf("Geo     : Opening geometry file: %s\n", geoFilePath);
     gmshOpen(geoFilePath, &ierr);

     if (ierr != 0) {
         fprintf(stderr, "Gmsh Error: Could not open GEO file '%s'. Code: %d\n", geoFilePath, ierr);
         ErrorGmsh(ierr);
         return; 
     }
 
     // Assigne la fonction de callback pour la taille de maillage
     geoSetSizeCallback(geoSize);
 
     printf("Geo     : Synchronizing OCC model...\n");
     gmshModelOccSynchronize(&ierr); ErrorGmsh(ierr);
 
     // Option Gmsh pour assurer la sauvegarde de toutes les entités (nœuds, etc.)
     gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr); ErrorGmsh(ierr);
 
     // Génération effective du maillage 2D
     printf("Geo     : Starting 2D mesh generation (using geoSize callback)...\n");
     gmshModelMeshGenerate(2, &ierr); ErrorGmsh(ierr);
     printf("Geo     : Mesh generated internally by Gmsh.\n");
 }
 
 



 // =========================================================================
 //        Section 2: Elasticity Problem Assembly and Solution
 // =========================================================================
 
 // Variables globales pour stocker une copie de A et B avant application des C.L. Dirichlet.
 // Nécessaires pour le calcul correct des forces résiduelles F_res = A*U - B.
 // Initialisées et libérées dans femElasticitySolve/Forces.
 static double **A_copy = NULL;
 static double *B_copy  = NULL;
 
 /**
  * @brief Assemble la contribution des éléments 2D à la matrice de rigidité (A) et au vecteur force (B).
  *
  * Itère sur tous les éléments 2D du maillage. Pour chaque élément :
  * - Itère sur les points d'intégration de Gauss.
  * - Calcule le Jacobien de la transformation et les dérivées des fonctions de forme.
  * - Calcule les contributions élémentaires à la matrice A (termes A, B, C de l'élasticité).
  * - Calcule la contribution élémentaire au vecteur B (terme de gravité).
  * - Ajoute ces contributions aux positions correspondantes dans le système global A et B.
  *
  * @param theProblem Pointeur vers la structure du problème contenant toutes les informations nécessaires.
  */
 void femElasticityAssembleElements(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;

    double x[3], y[3], phi[3], dphidxsi[3], dphideta[3], dphidx[3], dphidy[3];
    int iElem, iInteg, i, j, map[3], mapX[3], mapY[3];

    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j]; //{3, 8, 4}, map[0]=3, map[1]=8, map[2]=4.
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            // Use femDiscretePhi2 and femDiscreteDphi2
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
                }
                B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        }
    }         
}
 

 /**
  * @brief Assemble la contribution des conditions aux limites de Neumann au vecteur force (B).
  *
  * Itère sur toutes les conditions aux limites définies dans le problème.
  * Si une condition est de type NEUMANN_X ou NEUMANN_Y :
  * - Itère sur toutes les arêtes (éléments 1D) appartenant au domaine de cette condition.
  * - Pour chaque arête :
  *   - Itère sur les points d'intégration de Gauss sur l'arête.
  *   - Calcule le Jacobien de l'arête (demi-longueur).
  *   - Évalue les fonctions de forme 1D.
  *   - Ajoute la contribution de la force linéique imposée (value) au vecteur B global,
  *     aux DDL (Ux ou Uy) des nœuds de l'arête.
  *
  * @param theProblem Pointeur vers la structure du problème.
  */
 void femElasticityAssembleNeumann(femProblem *theProblem) {
     // Raccourcis
     femFullSystem  *theSystem = theProblem->system;
     femIntegration *theRule = theProblem->ruleEdge; // Règle d'intégration sur les arêtes
     femDiscrete    *theSpace = theProblem->spaceEdge;// Espace EF sur les arêtes (P1)
     femGeo         *theGeometry = theProblem->geometry;
     femNodes       *theNodes = theGeometry->theNodes;
     femMesh        *theEdges = theGeometry->theEdges; // Maillage contenant TOUTES les arêtes
 
     // Vérifications
     if (!theSystem || !theRule || !theSpace || !theGeometry || !theNodes ) {
         Error("AssembleNeumann: Missing required data structures."); return; }
     if (!theEdges) {
       printf("Warning: AssembleNeumann: No edges found in geometry, Neumann conditions will be ignored.\n"); return; }
     if (theEdges->nLocalNode != 2) {
         Error("AssembleNeumann: Edges mesh does not have 2 local nodes."); return; }
 
     double x[2], y[2], phi[2]; // Tableaux locaux pour une arête (2 nœuds)
     int map[2], mapU[2];       // Mapping local -> global
     int nLocal = 2;
     double *B = theSystem->B;
 
     // Boucle sur toutes les conditions aux limites définies
     for (int iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
         femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
         if (!theCondition) continue; // Skip si condition invalide
         femBoundaryType type = theCondition->type;
 
         // Sélectionne uniquement les conditions de Neumann
         int shift = -1; // 0 pour X, 1 pour Y
         if (type == NEUMANN_X) shift = 0;
         else if (type == NEUMANN_Y) shift = 1;
         else continue; // Passe à la condition suivante si ce n'est pas Neumann
 
         double value = theCondition->value;           // Valeur de la force linéique [N/m]
         femDomain *theDomain = theCondition->domain; // Domaine géométrique où s'applique la CL
         if (!theDomain || !theDomain->elem || theDomain->nElem == 0) continue; // Skip si domaine invalide ou vide
 
         // Boucle sur les arêtes appartenant à ce domaine
         for (int iEdgeInDomain = 0; iEdgeInDomain < theDomain->nElem; iEdgeInDomain++) {
             int iEdgeGlobal = theDomain->elem[iEdgeInDomain]; // Index de l'arête dans theEdges->elem
 
             // Vérification de la validité de l'index global de l'arête
             if (iEdgeGlobal < 0 || iEdgeGlobal >= theEdges->nElem) {
                  fprintf(stderr,"Warning: AssembleNeumann: Invalid edge index %d found in domain '%s'.\n", iEdgeGlobal, theDomain->name);
                  continue; // Passe à l'arête suivante dans ce domaine
             }
 
             // Récupération des coordonnées et mapping des DDL pour l'arête courante
             bool edgeValid = true;
             for (int j = 0; j < nLocal; j++) {
                 map[j] = theEdges->elem[iEdgeGlobal * nLocal + j]; // Indice global du nœud j de l'arête
                 if (map[j] == -1) {
                    edgeValid = false;
                    break; // Sortir de la boucle interne (j)
                }
                 // Vérification cruciale: le nœud existe-t-il encore après geoMeshFix?
                 if (map[j] < 0 || map[j] >= theNodes->nNodes) {
                      fprintf(stderr,"Warning: AssembleNeumann: Edge %d in domain '%s' uses invalid node index %d (possibly removed by geoMeshFix).\n", iEdgeGlobal, theDomain->name, map[j]);
                      edgeValid = false; break; // Ne pas traiter cette arête
                 }
                 mapU[j] = 2 * map[j] + shift; // Indice global du DDL (Ux ou Uy)
                 x[j] = theNodes->X[map[j]];
                 y[j] = theNodes->Y[map[j]];
             }
             if (!edgeValid) continue; // Passe à l'arête suivante si un nœud est invalide
 
             // Calcul du Jacobien de l'arête (demi-longueur)
             double edgeLength = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2));
             double jacEdge = edgeLength / 2.0;
             if (jacEdge < 1e-12) { // Vérifie si l'arête a une longueur nulle
                  fprintf(stderr,"Warning: AssembleNeumann: Zero length edge %d in domain '%s'.\n", iEdgeGlobal, theDomain->name);
                  continue; // Passe à l'arête suivante
             }
 
             // Intégration de Gauss sur l'arête
             for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
                 double xsi = theRule->xsi[iInteg];     // Point de Gauss sur [-1, 1]
                 double weight = theRule->weight[iInteg]; // Poids de Gauss
 
                 // Évaluation des fonctions de forme 1D au point de Gauss
                 femDiscretePhi(theSpace, xsi, phi);
 
                 // Ajout de la contribution au vecteur B global pour chaque nœud de l'arête
                 for (int i = 0; i < nLocal; i++) {
                     B[mapU[i]] += jacEdge * weight * phi[i] * value;
                 }
             } // Fin boucle points d'intégration
         } // Fin boucle arêtes du domaine
     } // Fin boucle conditions aux limites
 }
 
 
 /**
  * @brief Résout le problème d'élasticité linéaire.
  *
  * Orchestre les étapes de la résolution :
  * 1. Initialise le système linéaire (met A et B à zéro).
  * 2. Assemble les contributions des éléments 2D (matrice A, forces volumiques dans B).
  * 3. Assemble les contributions des conditions de Neumann (forces surfaciques/linéiques dans B).
  * 4. Sauvegarde une copie de A et B (avant application des C.L. Dirichlet) pour le calcul ultérieur des résidus.
  * 5. Applique les conditions de Dirichlet en modifiant A et B.
  * 6. Résout le système linéaire A*U = B par élimination de Gauss (via `femFullSystemEliminate`).
  * 7. Copie la solution U (qui se trouve dans B après élimination) dans `theProblem->soluce`.
  *
  * @param theProblem Pointeur vers la structure du problème.
  * @return double* Pointeur vers le vecteur solution (`theProblem->soluce`) contenant les déplacements nodaux (U).
  */
 double *femElasticitySolve(femProblem *theProblem) {
     femFullSystem *theSystem = theProblem->system;
     if (!theProblem || !theSystem) { Error("femElasticitySolve: Problem or system is NULL."); return NULL; }
 
     // 1. Initialisation (remise à zéro de A et B)
     femFullSystemInit(theSystem);
 
     // 2. Assemblage éléments (calcul de A, contribution gravité à B)
     femElasticityAssembleElements(theProblem);
 
     // 3. Assemblage Neumann (contribution forces surfaciques/linéiques à B)
     femElasticityAssembleNeumann(theProblem);
 
     int size = theSystem->size;
 
     // 4. Copie de A et B avant application de Dirichlet (pour calcul des résidus F = KU - B)
      if (size > 0) {
         // Allocation des copies si elles n'existent pas déjà
         if (A_copy == NULL) {
             A_copy = malloc(size * sizeof(double *));
             if (!A_copy) Error("Memory allocation failed for A_copy rows");
             // Allocation contiguë pour les données de A_copy
             A_copy[0] = malloc(size * size * sizeof(double));
              if (!A_copy[0]) { free(A_copy); A_copy=NULL; Error("Memory allocation failed for A_copy data"); }
             for (int i = 1; i < size; i++) A_copy[i] = A_copy[i-1] + size;
         }
         if (B_copy == NULL) {
             B_copy = malloc(size * sizeof(double));
             if (!B_copy) Error("Memory allocation failed for B_copy");
         }
         // Copie effective des valeurs actuelles de A et B
          for(int i=0; i<size; ++i)
             for(int j=0; j<size; ++j)
                 A_copy[i][j] = theSystem->A[i][j]; // Copie élément par élément
         memcpy(B_copy, theSystem->B, size * sizeof(double)); // Copie du vecteur B
      } else {
          // Si le système est vide, s'assurer que les copies sont NULL
          A_copy = NULL;
          B_copy = NULL;
      }
 
     // 5. Application des conditions de Dirichlet (modifie A et B dans theSystem)
     int *theConstrainedNodes = theProblem->constrainedNodes;
     if (theConstrainedNodes) {
         for (int i = 0; i < size; i++) {
             if (theConstrainedNodes[i] != -1) { // Si le DDL i est contraint
                 int conditionIndex = theConstrainedNodes[i];
                  // Vérifie la validité de l'index de condition
                  if (conditionIndex >= 0 && conditionIndex < theProblem->nBoundaryConditions) {
                     double value = theProblem->conditions[conditionIndex]->value;
                     femFullSystemConstrain(theSystem, i, value); // Applique la contrainte
                  } else {
                      fprintf(stderr,"Warning: Invalid condition index %d for constrained DOF %d\n", conditionIndex, i);
                  }
             }
         }
     }
 
     // 6. Résolution du système linéaire modifié A*U = B par élimination de Gauss
     //    La solution U écrase B dans theSystem->B.
     femFullSystemEliminate(theSystem);
 
     // 7. Copie de la solution U (qui est maintenant dans theSystem->B) vers theProblem->soluce
     if (size > 0 && theProblem->soluce) {
       memcpy(theProblem->soluce, theSystem->B, size * sizeof(double));
     }
 
     // Retourne un pointeur vers le vecteur solution stocké dans le problème
     return theProblem->soluce;
 }
 
 /**
  * @brief Calcule les forces résiduelles nodales F_res = A*U - B.
  *
  * Utilise la solution U calculée (`theProblem->soluce`) et les copies
  * de la matrice A et du vecteur B *avant* l'application des conditions de Dirichlet
  * (stockées dans `A_copy` et `B_copy`). Le résultat est stocké dans `theProblem->residuals`.
  * Les copies `A_copy` et `B_copy` sont ensuite libérées.
  *
  * @param theProblem Pointeur vers la structure du problème.
  * @return double* Pointeur vers le vecteur des forces résiduelles (`theProblem->residuals`).
  */
 double *femElasticityForces(femProblem *theProblem) {
     if (!theProblem) { Error("femElasticityForces: Problem is NULL."); return NULL; }
 
     double *residuals = theProblem->residuals;
     double *soluce = theProblem->soluce;
     int size = theProblem->system->size;
 
     // Alloue le vecteur residuals s'il n'existe pas (ne devrait pas arriver si Create est bien fait)
     if (residuals == NULL && size > 0) {
         residuals = malloc(size * sizeof(double));
         if (!residuals) Error("Failed to allocate residuals vector in femElasticityForces");
         theProblem->residuals = residuals;
     }
 
     // Si le système est vide, retourne le vecteur (vide ou NULL)
     if (size <= 0) return residuals;
 
     // Vérifie que les copies A_copy et B_copy existent (elles ont dû être créées dans Solve)
     if (A_copy == NULL || B_copy == NULL) {
         fprintf(stderr, "Warning: A_copy or B_copy is NULL in femElasticityForces. Cannot compute residuals accurately. Residuals set to zero.\n");
          if(residuals) { // Met à zéro si le vecteur existe
              for (int i = 0; i < size; i++) residuals[i] = 0.0;
          }
          // Tente de libérer au cas où l'un existe et pas l'autre (situation anormale)
          if (A_copy) { if(A_copy[0]) free(A_copy[0]); free(A_copy); A_copy = NULL; }
          if (B_copy) { free(B_copy); B_copy = NULL; }
         return residuals;
     }
 
     // Vérifie que soluce et residuals sont valides
      if (!soluce || !residuals) {
          Error("femElasticityForces: soluce or residuals vector is NULL.");
          // Libérer A_copy et B_copy avant de retourner
           if (A_copy) { if(A_copy[0]) free(A_copy[0]); free(A_copy); A_copy = NULL; }
           if (B_copy) { free(B_copy); B_copy = NULL; }
          return residuals; // Ou NULL ?
      }
 
     // Calcul effectif des résidus : R = A_copy * U - B_copy
     for (int i = 0; i < size; i++) {
         residuals[i] = -B_copy[i]; // Commence par -B[i]
         for (int j = 0; j < size; j++) {
             residuals[i] += A_copy[i][j] * soluce[j]; // Ajoute A[i][j] * U[j]
         }
     }
 
     // Libération impérative de la mémoire des copies A_copy et B_copy
     // pour éviter les fuites si Solve/Forces est appelé plusieurs fois.
     if (A_copy != NULL) {
         if (A_copy[0] != NULL) free(A_copy[0]); // Libère le bloc de données contigu
         free(A_copy);                           // Libère le tableau de pointeurs
         A_copy = NULL;                          // Réinitialise le pointeur global
     }
     if (B_copy != NULL) {
         free(B_copy);
         B_copy = NULL;
     }
 
     // Retourne le pointeur vers le vecteur de résidus stocké dans le problème
     return residuals;
 }
 



 
 // =========================================================================
 //                 Section 3: Utility Functions (Debugging)
 // =========================================================================
 
 /**
  * @brief Vérifie si une matrice carrée est symétrique à une tolérance près.
  *
  * @param matrix Matrice à vérifier (double**).
  * @param size Taille de la matrice (nombre de lignes/colonnes).
  * @param epsilon Tolérance absolue pour la comparaison (fabs(A[i][j] - A[j][i]) <= epsilon).
  * @return int 1 si symétrique, 0 sinon.
  */
 int isSymmetrical(double **matrix, int size, double epsilon) {
     if (!matrix) return 0; // Considère une matrice NULL comme non symétrique
     for (int i = 0; i < size; i++) {
         // Vérifie si la ligne existe (sécurité)
         if (!matrix[i]) return 0;
         // Compare A[i][j] et A[j][i] pour j > i
         for (int j = i + 1; j < size; j++) {
              // Vérifie si la ligne j existe (sécurité)
              if (!matrix[j]) return 0;
             if (fabs(matrix[i][j] - matrix[j][i]) > epsilon) {
                 // Optionnel: Afficher où la symétrie échoue
                 // printf("Non-symmetry detected: A[%d][%d] = %e != A[%d][%d] = %e\n", i, j, matrix[i][j], j, i, matrix[j][i]);
                 return 0; // Non symétrique
             }
         }
     }
     return 1; // Symétrique
 }
 
 /**
  * @brief Affiche une matrice carrée de manière formatée (en cachant les zéros proches).
  *
  * Affiche seulement les `maxDisplay` premières lignes/colonnes si la matrice est grande.
  * Remplace les valeurs proches de zéro (selon `epsilon`) par des espaces pour la lisibilité.
  *
  * @param A Matrice à afficher.
  * @param size Taille de la matrice.
  * @param maxDisplay Nombre maximum de lignes/colonnes à afficher.
  * @param epsilon Seuil pour considérer une valeur comme nulle.
  */
 void printMatrixClean(double **A, int size, int maxDisplay, double epsilon) {
     printf("Affichage Matrice A [%d x %d] (max %d x %d, tol %.1e)\n", size, size, maxDisplay, maxDisplay, epsilon);
     printf("======================================================\n");
     if (!A || size <= 0) { printf("  (Matrice vide ou NULL)\n"); printf("======================================================\n\n"); return; }
 
     int limit = (size < maxDisplay) ? size : maxDisplay;
 
     for (int i = 0; i < limit; i++) {
         if (!A[i]) { printf("  (Ligne %d est NULL)\n", i); continue; } // Sécurité
         printf("  "); // Petite marge
         for (int j = 0; j < limit; j++) {
             if (fabs(A[i][j]) > epsilon) {
                 printf("%+11.4e ", A[i][j]); // Format scientifique aligné
             } else {
                 printf("            ");     // 12 espaces pour l'alignement
             }
         }
         // Afficher "..." si toutes les colonnes ne sont pas montrées
         if (limit < size) printf(" ...");
         printf("\n");
     }
     // Afficher "..." si toutes les lignes ne sont pas montrées
     if (limit < size) {
         printf("  ...\n");
         printf("  (Seulement les %d premières lignes/colonnes affichées)\n", maxDisplay);
     }
     printf("======================================================\n\n");
 }
 
 /**
  * @brief Affiche les éléments diagonaux d'une matrice carrée.
  *
  * @param A Matrice dont afficher la diagonale.
  * @param size Taille de la matrice.
  */
 void printDiagonal(double **A, int size) {
     printf("Affichage Diagonale Matrice A [%d x %d]\n", size, size);
     printf("======================================================\n");
      if (!A || size <= 0) { printf("  (Matrice vide ou NULL)\n"); printf("======================================================\n\n"); return; }
 
     for (int i = 0; i < size; i++) {
          if (!A[i]) { printf("  A[%d][%d] : Ligne NULL\n", i, i); continue;} // Sécurité
         printf("  A[%d][%d] = %12.4e\n", i, i, A[i][i]);
     }
     printf("======================================================\n\n");
 }