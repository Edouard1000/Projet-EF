#include "fem.h"
#include <math.h>
#include <stdbool.h>


// === Fonction pour l'interpolation de la taille du maillage (Mesh) ===

double hermiteInterpolation(double d, double d_star, double h0, double h_star) {
    if (d >= d_star) return h_star;
    double t = d / d_star;
    return h0 + (h_star - h0) * (3 * t * t - 2 * t * t * t);
}

double geoSize(double x, double y){
    femGeo* theGeometry = geoGetGeometry();
    double h = theGeometry->h;

    // Trouver le sommet du rail et points intermédiaires
    static double yTop, X1, Y1, X2, Y2;
    yTop = 172;
    X1 = 37.85;
    Y1 = 157.7;
    X2 = -X1;
    Y2 = Y1;

    // Calcul des distances
    double dTop = fabs(yTop - y);
    double d1 = sqrt((x - X1) * (x - X1) + (y - Y1) * (y - Y1));
    double d2 = sqrt((x - X2) * (x - X2) + (y - Y2) * (y - Y2));

    // Tailles de mailles plus fines
    double hTop = h * 0.1;
    double hUpper = hermiteInterpolation(dTop, 35, hTop, h);
    double hPoint1 = hermiteInterpolation(d1, 25, hTop, h);
    double hPoint2 = hermiteInterpolation(d2, 25, hTop, h);

    return fmin(fmin(hUpper, hPoint1), hPoint2); // Retourne la plus petite taille
}


// === Assemblage des éléments (LinearElasticity) ===


// Function to check if a matrix is symmetrical
bool isSymmetrical(double **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

void femElasticityAssembleElements(femProblem *theProblem) {
  // Récupération des structures de données
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,i,j,map[4],mapX[4],mapY[4];

    // Paramètres du problème
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
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        }

        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i]*dphidxsi[i];
                dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i];
                dydeta += y[i]*dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                                            dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                                            dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                                            dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                                            dphidx[i] * c * dphidx[j]) * jac * weight;
                }
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; // Ajout du terme source
            }
        }
    }

    // Check if the matrix is symmetrical
    if (isSymmetrical(theProblem->system->A, theProblem->system->size)) {
        printf("The matrix is symmetrical.\n");
    } else {
        printf("The matrix is not symmetrical.\n");
    }
         
}


// Variables globales pour stocker les copies de A et B (pour le calcul des résidus)
double **A_copy = NULL;
double *B_copy  = NULL;

void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->ruleEdge;
    femDiscrete    *theSpace    = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theEdges    = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, map[2], mapU[2];

    int nLocal = 2;
    double *B  = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value = theCondition->value;

        if (type == DIRICHLET_X || type == DIRICHLET_Y) { continue; }   // Ignore les Dirichlet

        int shift = (type == NEUMANN_X) ? 0 : 1;                        // Décalage pour X ou Y

        for (iEdge = 0; iEdge < domain->nElem; iEdge++)
        {
            iElem = domain->elem[iEdge];
            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }

            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx * dx + dy * dy);
            double jac = length / 2;

            for (iInteg = 0; iInteg < theRule->n; iInteg++)
            {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];
                femDiscretePhi(theSpace, xsi, phi);
                for (i = 0; i < theSpace->n; i++) {
                    B[mapU[i]] += phi[i] * value * jac * weight;
                }
            }
        }
    }
}


double *femElasticitySolve(femProblem *theProblem) {
    femFullSystem *theSystem = theProblem->system;
    femFullSystemInit(theSystem);

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    int size = theSystem->size;

    // Allocation et copie de A et B pour le calcul des résidus
    if (A_copy == NULL) {
        A_copy = (double **)malloc(sizeof(double *) * size);
        for (int i = 0; i < size; i++) {
            A_copy[i] = (double *)malloc(sizeof(double) * size);
        }
    }
    if (B_copy == NULL) {
        B_copy = (double *)malloc(sizeof(double) * size);
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            A_copy[i][j] = theSystem->A[i][j];
        }
        B_copy[i] = theSystem->B[i];
    }

    // Application des conditions de Dirichlet
    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    // Résolution
    femFullSystemEliminate(theSystem);
    memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    return theProblem->soluce;
}


double *femElasticityForces(femProblem *theProblem) {
    double *residuals = theProblem->residuals;
    double *soluce = theProblem->soluce;
    int size = theProblem->system->size;

    if (residuals == NULL) {
        residuals = (double *)malloc(sizeof(double) * size);
    }
    for (int i = 0; i < size; i++) {
        residuals[i] = 0.0;
    }

    // Calcul des résidus : R = A * U - B
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            residuals[i] += A_copy[i][j] * soluce[j];
        }
        residuals[i] -= B_copy[i];
    }

    // Libération de la mémoire des copies (IMPORTANT : évite les fuites de mémoire)
    for (int i = 0; i < size; i++) {
      free(A_copy[i]);
      A_copy[i] = NULL;
    }

    free(A_copy);
    free(B_copy);
    A_copy = NULL;
    B_copy = NULL;

    return residuals;
}



void geoMeshGenerate() {
    int ierr;

    femGeo* theGeometry = geoGetGeometry();
    
    gmshOpen("../Projet-EF/UIC60_copie2.geo", &ierr); //OUVRE LE FICHIER .GEO
    ErrorGmsh(ierr);

    geoSetSizeCallback(geoSize);

    gmshModelOccSynchronize(&ierr); //SYNCHRONISE LE MODÈLE

    // Option pour que Gmsh sauvegarde *tous* les éléments (pas seulement ceux avec des propriétés physiques)
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);

     // Générer le maillage 2D
    gmshModelMeshGenerate(2, &ierr);
}