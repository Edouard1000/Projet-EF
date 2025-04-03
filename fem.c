/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

 #include "fem.h"
 #include <stdbool.h>


 femGeo theGeometry;
 
 femGeo *geoGetGeometry()                        { return &theGeometry; }
 
 double geoSizeDefault(double x, double y)       { return theGeometry.h; }
 
 double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data)
                                                 { return theGeometry.geoSize(x,y);    }
 void geoInitialize() 
 {
     int ierr;
     theGeometry.geoSize = geoSizeDefault;
     gmshInitialize(0,NULL,1,0,&ierr);                         ErrorGmsh(ierr);
     gmshModelAdd("MyGeometry",&ierr);                         ErrorGmsh(ierr);
     gmshModelMeshSetSizeCallback(geoGmshSize,NULL,&ierr);     ErrorGmsh(ierr);
     theGeometry.theNodes = NULL;
     theGeometry.theElements = NULL;
     theGeometry.theEdges = NULL;
     theGeometry.nDomains = 0;
     theGeometry.theDomains = NULL;
 }
 
 void geoFinalize() 
 {
     int ierr;
     
     if (theGeometry.theNodes) {
         free(theGeometry.theNodes->X);
         free(theGeometry.theNodes->Y);
         free(theGeometry.theNodes); }
     if (theGeometry.theElements) {
         free(theGeometry.theElements->elem);
         free(theGeometry.theElements); }
     if (theGeometry.theEdges) {
         free(theGeometry.theEdges->elem);
         free(theGeometry.theEdges); }
     for (int i=0; i < theGeometry.nDomains; i++) {
         free(theGeometry.theDomains[i]->elem);
         free(theGeometry.theDomains[i]);  }
     free(theGeometry.theDomains);
     gmshFinalize(&ierr); ErrorGmsh(ierr);
 }
 
 
 void geoSetSizeCallback(double (*geoSize)(double x, double y)) 
 {
     theGeometry.geoSize = geoSize; }
 
 
 void geoMeshImport() 
 {
     int ierr;
     
     /* Importing nodes */
     
     size_t nNode,n,m,*node;
     double *xyz,*trash;
     gmshModelMeshGetNodes(&node,&nNode,&xyz,&n,
                          &trash,&m,-1,-1,0,0,&ierr);          ErrorGmsh(ierr);                         
     femNodes *theNodes = malloc(sizeof(femNodes));
     theNodes->nNodes = nNode;
     theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
     theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
     for (int i = 0; i < theNodes->nNodes; i++){
         theNodes->X[i] = xyz[3*node[i]-3];
         theNodes->Y[i] = xyz[3*node[i]-2]; }
     theGeometry.theNodes = theNodes;
     gmshFree(node);
     gmshFree(xyz);
     gmshFree(trash);
     printf("Geo     : Importing %d nodes \n",theGeometry.theNodes->nNodes);
        
     /* Importing elements */
     /* Pas super joli : a ameliorer pour eviter la triple copie */
         
     size_t nElem, *elem;
     gmshModelMeshGetElementsByType(1,&elem,&nElem,
                                &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
     femMesh *theEdges = malloc(sizeof(femMesh));
     theEdges->nLocalNode = 2;
     theEdges->nodes = theNodes;
     theEdges->nElem = nElem;  
     theEdges->elem = malloc(sizeof(int)*2*theEdges->nElem);
     for (int i = 0; i < theEdges->nElem; i++)
         for (int j = 0; j < theEdges->nLocalNode; j++)
             theEdges->elem[2*i+j] = node[2*i+j]-1;  
     theGeometry.theEdges = theEdges;
     int shiftEdges = elem[0];
     gmshFree(node);
     gmshFree(elem);
     printf("Geo     : Importing %d edges \n",theEdges->nElem);
   
     gmshModelMeshGetElementsByType(2,&elem,&nElem,
                                &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
     if (nElem != 0) {
       femMesh *theElements = malloc(sizeof(femMesh));
       theElements->nLocalNode = 3;
       theElements->nodes = theNodes;
       theElements->nElem = nElem;  
       theElements->elem = malloc(sizeof(int)*3*theElements->nElem);
       for (int i = 0; i < theElements->nElem; i++)
           for (int j = 0; j < theElements->nLocalNode; j++)
               theElements->elem[3*i+j] = node[3*i+j]-1;  
       theGeometry.theElements = theElements;
       gmshFree(node);
       gmshFree(elem);
       printf("Geo     : Importing %d triangles \n",theElements->nElem); }
     
     int nElemTriangles = nElem;
     gmshModelMeshGetElementsByType(3,&elem,&nElem,
                                &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
     if (nElem != 0 && nElemTriangles != 0)  
       Error("Cannot consider hybrid geometry with triangles and quads :-(");                       
                                
     if (nElem != 0) {
       femMesh *theElements = malloc(sizeof(femMesh));
       theElements->nLocalNode = 4;
       theElements->nodes = theNodes;
       theElements->nElem = nElem;  
       theElements->elem = malloc(sizeof(int)*4*theElements->nElem);
       for (int i = 0; i < theElements->nElem; i++)
           for (int j = 0; j < theElements->nLocalNode; j++)
               theElements->elem[4*i+j] = node[4*i+j]-1;  
       theGeometry.theElements = theElements;
       gmshFree(node);
       gmshFree(elem);
       printf("Geo     : Importing %d quads \n",theElements->nElem); }
 
     
     /* Importing 1D entities */
   
     int *dimTags;
     gmshModelGetEntities(&dimTags,&n,1,&ierr);        ErrorGmsh(ierr);
     theGeometry.nDomains = n/2;
     theGeometry.theDomains = malloc(sizeof(femDomain*)*n/2);
     printf("Geo     : Importing %d entities \n",theGeometry.nDomains);
 
     for (int i=0; i < n/2; i++) {
         int dim = dimTags[2*i+0];
         int tag = dimTags[2*i+1];
         femDomain *theDomain = malloc(sizeof(femDomain)); 
         theGeometry.theDomains[i] = theDomain;
         theDomain->mesh = theEdges;
         sprintf(theDomain->name, "Entity %d ",tag-1);
          
         int *elementType;
         size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags; 
         gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
         theDomain->nElem = nElementTags[0];
         theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
         for (int j = 0; j < theDomain->nElem; j++) {
             theDomain->elem[j] = elementTags[0][j] - shiftEdges; }
         printf("Geo     : Entity %d : %d elements \n",i,theDomain->nElem);
         gmshFree(nElementTags);
         gmshFree(nNodesTags);
         gmshFree(elementTags);
         gmshFree(nodesTags);
         gmshFree(elementType); }
     gmshFree(dimTags);
  
     return;
 
 }
 
 void geoMeshPrint() 
 {
    femNodes *theNodes = theGeometry.theNodes;
    if (theNodes != NULL) {
       printf("Number of nodes %d \n", theNodes->nNodes);
       for (int i = 0; i < theNodes->nNodes; i++) {
         printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }}
    femMesh *theEdges = theGeometry.theEdges;
    if (theEdges != NULL) {
      printf("Number of edges %d \n", theEdges->nElem);
      int *elem = theEdges->elem;
      for (int i = 0; i < theEdges->nElem; i++) {
         printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
    femMesh *theElements = theGeometry.theElements;
    if (theElements != NULL) {
      if (theElements->nLocalNode == 3) {
         printf("Number of triangles %d \n", theElements->nElem);
         int *elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
      if (theElements->nLocalNode == 4) {
         printf("Number of quads %d \n", theElements->nElem);
         int *elem = theElements->elem;
         for (int i = 0; i < theElements->nElem; i++) {
             printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
    int nDomains = theGeometry.nDomains;
    printf("Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = theGeometry.theDomains[iDomain];
       printf("  Domain : %6d \n", iDomain);
       printf("  Name : %s\n", theDomain->name);
       printf("  Number of elements : %6d\n", theDomain->nElem);
       for (int i=0; i < theDomain->nElem; i++){
  //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
           printf("%6d",theDomain->elem[i]);
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
       printf("\n"); }
   
   
 }
 
 
 void geoMeshWrite(const char *filename) 
 {
    FILE* file = fopen(filename,"w");
  
    femNodes *theNodes = theGeometry.theNodes;
    fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) {
       fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
       
    femMesh *theEdges = theGeometry.theEdges;
    fprintf(file,"Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++) {
       fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
       
    femMesh *theElements = theGeometry.theElements;
    if (theElements->nLocalNode == 3) {
       fprintf(file,"Number of triangles %d \n", theElements->nElem);
       elem = theElements->elem;
       for (int i = 0; i < theElements->nElem; i++) {
           fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
    if (theElements->nLocalNode == 4) {
       fprintf(file,"Number of quads %d \n", theElements->nElem);
       elem = theElements->elem;
       for (int i = 0; i < theElements->nElem; i++) {
           fprintf(file,"%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}
      
    int nDomains = theGeometry.nDomains;
    fprintf(file,"Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = theGeometry.theDomains[iDomain];
       fprintf(file,"  Domain : %6d \n", iDomain);
       fprintf(file,"  Name : %s\n", theDomain->name);
       fprintf(file,"  Number of elements : %6d\n", theDomain->nElem);
       for (int i=0; i < theDomain->nElem; i++){
           fprintf(file,"%6d",theDomain->elem[i]);
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
       fprintf(file,"\n"); }
     
    fclose(file);
 }
 
 void geoMeshRead(const char *filename) 
 {
    FILE* file = fopen(filename,"r");
    
    int trash, *elem;
    
    femNodes *theNodes = malloc(sizeof(femNodes));
    theGeometry.theNodes = theNodes;
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 
 
    femMesh *theEdges = malloc(sizeof(femMesh));
    theGeometry.theEdges = theEdges;
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
    theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
    for(int i=0; i < theEdges->nElem; ++i) {
         elem = theEdges->elem;
         ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
   
    femMesh *theElements = malloc(sizeof(femMesh));
    theGeometry.theElements = theElements;
    theElements->nLocalNode = 0;
    theElements->nodes = theNodes;
    char elementType[MAXNAME];  
    ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
    if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
       theElements->nLocalNode = 3;
       theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
       for(int i=0; i < theElements->nElem; ++i) {
           elem = theElements->elem;
           ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                     &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
    if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
       theElements->nLocalNode = 4;
       theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
       for(int i=0; i < theElements->nElem; ++i) {
           elem = theElements->elem;
           ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                     &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
            
    ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
    int nDomains = theGeometry.nDomains;
    theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
       femDomain *theDomain = malloc(sizeof(femDomain)); 
       theGeometry.theDomains[iDomain] = theDomain;
       theDomain->mesh = theEdges; 
       ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
       ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
       ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
       theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
       for (int i=0; i < theDomain->nElem; i++){
           ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
           if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
     
    fclose(file);
 }


 // =========================================================================
 //              Mesh Correction Function (Student Project Addition)
 // =========================================================================

 /**
  * @brief Implémentation de la correction du maillage (suppression des nœuds isolés).
  * @cite maillage.py by Vincent Legat, UCL-IMMC
  *
  * Fonction ajoutée pour reproduire le comportement d'un script externe (maillage.py)
  * qui supprimait les nœuds non référencés par les éléments 2D, corrigeant
  * potentiellement des problèmes de pivots nuls ou de maillage invalide généré
  * par certaines versions/configurations de Gmsh.
  *
  * Étapes :
  * 1. Parcours des éléments 2D pour marquer tous les nœuds utilisés.
  * 2. Comptage des nœuds utilisés et création d'une table de renumérotation
  *    (ancien indice -> nouvel indice, ou -1 si supprimé).
  * 3. Création de nouveaux tableaux pour les coordonnées X et Y contenant
  *    uniquement les nœuds conservés.
  * 4. Création d'un nouveau tableau pour la connectivité des éléments 2D
  *    avec les indices de nœuds mis à jour selon la table de renumérotation.
  * 5. Remplacement des anciens tableaux (coordonnées, connectivité 2D) par les
  *    nouveaux dans la structure femGeo et mise à jour du nombre de nœuds.
  * 6. Libération de la mémoire temporaire utilisée.
  *
  * @param geometry Pointeur vers la structure femGeo à corriger en place.
  */

 void geoMeshFix(femGeo *geometry) {
    printf("Geo     : Running mesh fix...\n");
    if (!geometry || !geometry->theNodes || !geometry->theElements) {
        Warning("geoMeshFix called with invalid geometry data.");
        return;
    }

    femNodes *theNodes = geometry->theNodes;
    femMesh *theElements = geometry->theElements;
    // On ne touche PAS aux theEdges et theDomains ici pour préserver les indices

    int nnodes_old = theNodes->nNodes;
    int nelem = theElements->nElem;
    int nlocal_elem = theElements->nLocalNode;

    if (nnodes_old == 0 || nelem == 0) {
        printf("Geo     : Mesh fix skipped (no nodes or elements).\n");
        return;
    }

    // 1. Marquer les nœuds utilisés UNIQUEMENT par les éléments 2D
    bool *node_used = calloc(nnodes_old, sizeof(bool)); // Initialise à false
    if (!node_used) Error("Memory allocation failed for node_used in geoMeshFix");

    for (int i = 0; i < nelem; ++i) {
        for (int j = 0; j < nlocal_elem; ++j) {
            int node_idx = theElements->elem[i * nlocal_elem + j];
            if (node_idx >= 0 && node_idx < nnodes_old) {
                node_used[node_idx] = true;
            } else {
                 fprintf(stderr, "Warning: Invalid node index %d found in element %d during mesh fix.\n", node_idx, i);
            }
        }
    }

    // 2. Compter les nœuds utilisés et créer la table de renumérotation
    int nnodes_new = 0;
    for (int i = 0; i < nnodes_old; ++i) {
        if (node_used[i]) {
            nnodes_new++;
        }
    }

    if (nnodes_new == nnodes_old) {
        printf("Geo     : Mesh fix not needed. All nodes are used by 2D elements.\n");
        free(node_used);
        return; // Pas besoin de corriger
    }
     if (nnodes_new == 0) {
         printf("Geo     : Mesh fix resulted in zero nodes. Aborting fix.\n");
         free(node_used);
         return; // Évite de créer un maillage vide
     }

    printf("Geo     : Fixing mesh: Removing %d nodes not used by 2D elements.\n", nnodes_old - nnodes_new);

    int *renumber = malloc(nnodes_old * sizeof(int));
    if (!renumber) { free(node_used); Error("Memory allocation failed for renumber map in geoMeshFix"); }

    int current_new_index = 0;
    for (int i = 0; i < nnodes_old; ++i) {
        if (node_used[i]) {
            renumber[i] = current_new_index++;
        } else {
            renumber[i] = -1; // Marquer comme non utilisé/supprimé
        }
    }

    // 3. Créer les nouvelles listes de coordonnées
    double *newX = malloc(nnodes_new * sizeof(double));
    double *newY = malloc(nnodes_new * sizeof(double));
    if (!newX || !newY) {
        free(node_used); free(renumber);
        Error("Memory allocation failed for new node coordinates in geoMeshFix");
    }

    current_new_index = 0;
    for (int i = 0; i < nnodes_old; ++i) {
        if (node_used[i]) { // ou renumber[i] != -1
            newX[current_new_index] = theNodes->X[i];
            newY[current_new_index] = theNodes->Y[i];
            current_new_index++;
        }
    }

    // 4. Créer la nouvelle liste d'éléments 2D avec les indices renumérotés
    int *newElem = malloc(nelem * nlocal_elem * sizeof(int));
    if (!newElem) {
        free(node_used); free(renumber); free(newX); free(newY);
        Error("Memory allocation failed for new element connectivity in geoMeshFix");
    }
    for (int i = 0; i < nelem; ++i) {
        for (int j = 0; j < nlocal_elem; ++j) {
            int old_node_idx = theElements->elem[i * nlocal_elem + j];
            if (old_node_idx >= 0 && old_node_idx < nnodes_old) {
                 int new_node_idx = renumber[old_node_idx];
                  if (new_node_idx == -1) {
                      // Cela ne devrait pas arriver car on a marqué les noeuds utilisés par ces éléments
                       fprintf(stderr,"Error: Inconsistency in geoMeshFix - node %d used by element %d marked for removal.\n", old_node_idx, i);
                       newElem[i * nlocal_elem + j] = -1; // Indique une erreur
                  } else {
                       newElem[i * nlocal_elem + j] = new_node_idx;
                  }
            } else {
                 fprintf(stderr,"Warning: Invalid old node index %d in element %d during renumbering.\n", old_node_idx, i);
                 newElem[i * nlocal_elem + j] = -1; // Indique une erreur
            }
        }
    }

     // 5. Mettre à jour la connectivité des ARÊTES (theEdges)
    if (geometry->theEdges && geometry->theEdges->elem) {
        printf("    Updating edge connectivity with renumbered nodes...\n");
        femMesh *theEdges = geometry->theEdges;
        int nEdges = theEdges->nElem;
        int nLocalEdge = theEdges->nLocalNode; // Devrait être 2
        int *edgesElem = theEdges->elem;
        int invalid_node_count_in_edges = 0;

        for (int i = 0; i < nEdges * nLocalEdge; ++i) { // Parcours linéaire du tableau
            int old_node_idx = edgesElem[i];
            if (old_node_idx >= 0 && old_node_idx < nnodes_old) {
                int new_node_idx = renumber[old_node_idx];
                edgesElem[i] = new_node_idx; // Remplace par le nouvel indice ou -1
                if (new_node_idx == -1) {
                    invalid_node_count_in_edges++;
                }
            } else {
                // L'indice était déjà invalide avant la correction ?
                fprintf(stderr, "    Warning: Found invalid node index %d in original edgesElem[%d]. Setting to -1.\n", old_node_idx, i);
                edgesElem[i] = -1;
                invalid_node_count_in_edges++;
            }
        }
        if (invalid_node_count_in_edges > 0) {
            printf("    Info: %d node entries in edges now point to removed nodes (index -1).\n", invalid_node_count_in_edges);
        }
    } else {
        printf("    Skipping edge connectivity update (no edges or elem data).\n");
    }

    // 6. Libérer les anciennes données de nœuds et d'éléments 2D et mettre à jour la structure geometry
    free(theNodes->X);
    free(theNodes->Y);
    free(theElements->elem);

    theNodes->X = newX;
    theNodes->Y = newY;
    theNodes->nNodes = nnodes_new;

    theElements->elem = newElem;
    // theElements->nElem et nLocalNode ne changent pas

    // 7. Libérer la mémoire temporaire
    free(node_used);
    free(renumber);

    printf("Geo     : Mesh fix completed. New node count: %d\n", nnodes_new);
}



 
 void geoSetDomainName(int iDomain, char *name) 
 {
     if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
     if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
     sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
 } 
 
 int geoGetDomain(char *name)
 {
     int theIndex = -1;
     int nDomains = theGeometry.nDomains;
     for (int iDomain = 0; iDomain < nDomains; iDomain++) {
         femDomain *theDomain = theGeometry.theDomains[iDomain];
         if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
             theIndex = iDomain;  }
     return theIndex;
             
 }
 
 
 static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
 static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
 static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
 static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
 static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
 static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
 static const double _gaussEdge2Xsi[2]    = { 0.577350269189626,-0.577350269189626};
 static const double _gaussEdge2Weight[2] = { 1.000000000000000, 1.000000000000000};
 
 
 
 femIntegration *femIntegrationCreate(int n, femElementType type)
 {
     femIntegration *theRule = malloc(sizeof(femIntegration));
     if (type == FEM_QUAD && n == 4) {
         theRule->n      = 4;
         theRule->xsi    = _gaussQuad4Xsi;
         theRule->eta    = _gaussQuad4Eta;
         theRule->weight = _gaussQuad4Weight; }
     else if (type == FEM_TRIANGLE && n == 3) {
         theRule->n      = 3;
         theRule->xsi    = _gaussTri3Xsi;
         theRule->eta    = _gaussTri3Eta;
         theRule->weight = _gaussTri3Weight; }
     else if (type == FEM_EDGE && n == 2) {
         theRule->n      = 2;
         theRule->xsi    = _gaussEdge2Xsi;
         theRule->eta    = NULL;
         theRule->weight = _gaussEdge2Weight; }
     else Error("Cannot create such an integration rule !");
     return theRule; 
 }
 
 void femIntegrationFree(femIntegration *theRule)
 {
     free(theRule);
 }
 
 void _q1c0_x(double *xsi, double *eta) 
 {
     xsi[0] =  1.0;  eta[0] =  1.0;
     xsi[1] = -1.0;  eta[1] =  1.0;
     xsi[2] = -1.0;  eta[2] = -1.0;
     xsi[3] =  1.0;  eta[3] = -1.0;
 }
 
 void _q1c0_phi(double xsi, double eta, double *phi)
 {
     phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
     phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
     phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
     phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
 }
 
 void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     dphidxsi[0] =   (1.0 + eta) / 4.0;  
     dphidxsi[1] = - (1.0 + eta) / 4.0;
     dphidxsi[2] = - (1.0 - eta) / 4.0;
     dphidxsi[3] =   (1.0 - eta) / 4.0;
     dphideta[0] =   (1.0 + xsi) / 4.0;  
     dphideta[1] =   (1.0 - xsi) / 4.0;
     dphideta[2] = - (1.0 - xsi) / 4.0;
     dphideta[3] = - (1.0 + xsi) / 4.0;
 
 }
 
 void _p1c0_x(double *xsi, double *eta) 
 {
     xsi[0] =  0.0;  eta[0] =  0.0;
     xsi[1] =  1.0;  eta[1] =  0.0;
     xsi[2] =  0.0;  eta[2] =  1.0;
 }
 
 void _p1c0_phi(double xsi, double eta, double *phi)
 {
     phi[0] = 1 - xsi - eta;  
     phi[1] = xsi;
     phi[2] = eta;
 }
 
 void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     dphidxsi[0] = -1.0;  
     dphidxsi[1] =  1.0;
     dphidxsi[2] =  0.0;
     dphideta[0] = -1.0;  
     dphideta[1] =  0.0;
     dphideta[2] =  1.0;
 }
 
 void _e1c0_x(double *xsi) 
 {
     xsi[0] = -1.0;  
     xsi[1] =  1.0;  
 }
 
 void _e1c0_phi(double xsi,  double *phi)
 {
     phi[0] = (1 - xsi) / 2.0;  
     phi[1] = (1 + xsi) / 2.0;
 }
 
 void _e1c0_dphidx(double xsi, double *dphidxsi)
 {
     dphidxsi[0] = -0.5;  
     dphidxsi[1] =  0.5;
 }
 
 
 
 femDiscrete *femDiscreteCreate(int n, femElementType type)
 {
     femDiscrete *theSpace = malloc(sizeof(femDiscrete));
     theSpace->type = type;  
     theSpace->n = 0;
     theSpace->x = NULL;
     theSpace->phi = NULL;   
     theSpace->dphidx = NULL;    
     theSpace->x2 = NULL;    
     theSpace->phi2 = NULL;
     theSpace->dphi2dx = NULL;
  
     if (type == FEM_QUAD && n == 4) {
         theSpace->n       = 4;
         theSpace->x2      = _q1c0_x;
         theSpace->phi2    = _q1c0_phi;
         theSpace->dphi2dx = _q1c0_dphidx; }
     else if (type == FEM_TRIANGLE && n == 3) {
         theSpace->n       = 3;
         theSpace->x2      = _p1c0_x;
         theSpace->phi2    = _p1c0_phi;
         theSpace->dphi2dx = _p1c0_dphidx; }
     else if (type == FEM_EDGE && n == 2) {
         theSpace->n       = 2;
         theSpace->x       = _e1c0_x;
         theSpace->phi     = _e1c0_phi;
         theSpace->dphidx  = _e1c0_dphidx; }
     else Error("Cannot create such a discrete space !");
     return theSpace; 
 }
 
 void femDiscreteFree(femDiscrete *theSpace)
 {
     free(theSpace);
 }
 
 void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
 {
     mySpace->x2(xsi,eta);
 }
 
 void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
 {
     mySpace->phi2(xsi,eta,phi);
 }
 
 void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
 {
     mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
 }
 
 void femDiscreteXsi(femDiscrete* mySpace, double *xsi)
 {
     mySpace->x(xsi);
 }
 
 void femDiscretePhi(femDiscrete* mySpace, double xsi, double *phi)
 {
     mySpace->phi(xsi,phi);
 }
 
 void femDiscreteDphi(femDiscrete* mySpace, double xsi, double *dphidxsi)
 {
     mySpace->dphidx(xsi,dphidxsi);
 }
 
 void femDiscretePrint(femDiscrete *mySpace)
 {
     int i,j;
     int n = mySpace->n;
     double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
 
     if (mySpace->type == FEM_EDGE) {
         femDiscreteXsi(mySpace,xsi);
         for (i=0; i < n; i++) {           
             femDiscretePhi(mySpace,xsi[i],phi);
             femDiscreteDphi(mySpace,xsi[i],dphidxsi);
             for (j=0; j < n; j++)  {
                 printf("(xsi=%+.1f) : ",xsi[i]);
                 printf(" phi(%d)=%+.1f",j,phi[j]);  
                 printf("   dphidxsi(%d)=%+.1f \n",j,dphidxsi[j]); }
             printf(" \n"); }}
     
     if (mySpace->type == FEM_QUAD || mySpace->type == FEM_TRIANGLE) {
         femDiscreteXsi2(mySpace, xsi, eta);
         for (i = 0; i < n; i++)  {    
             femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
             femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);
             for (j = 0; j < n; j++) {  
                 printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);  
                 printf(" phi(%d)=%+.1f", j, phi[j]);
                 printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
                 printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]); }
             printf(" \n"); }}   
 }
 
 femFullSystem *femFullSystemCreate(int size)
 {
     femFullSystem *theSystem = malloc(sizeof(femFullSystem));
     femFullSystemAlloc(theSystem, size);
     femFullSystemInit(theSystem);
 
     return theSystem; 
 }
 
 void femFullSystemFree(femFullSystem *theSystem)
 {
     free(theSystem->A);
     free(theSystem->B);
     free(theSystem);
 }
 
 void femFullSystemAlloc(femFullSystem *mySystem, int size)
 {
     int i;  
     double *elem = malloc(sizeof(double) * size * (size+1)); 
     mySystem->A = malloc(sizeof(double*) * size); 
     mySystem->B = elem;
     mySystem->A[0] = elem + size;  
     mySystem->size = size;
     for (i=1 ; i < size ; i++) 
         mySystem->A[i] = mySystem->A[i-1] + size;
 }
 
 void femFullSystemInit(femFullSystem *mySystem)
 {
     int i,size = mySystem->size;
     for (i=0 ; i < size*(size+1) ; i++) 
         mySystem->B[i] = 0;}
 
 
 void femFullSystemPrint(femFullSystem *mySystem)
 {
     double  **A, *B;
     int     i, j, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
     
     for (i=0; i < size; i++) {
         for (j=0; j < size; j++)
             if (A[i][j] == 0)  printf("         ");   
             else               printf(" %+.1e",A[i][j]);
         printf(" :  %+.1e \n",B[i]); }
 }

 
 
 double* femFullSystemEliminate(femFullSystem *mySystem)
 {
     double  **A, *B, factor;
     int     i, j, k, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
     
     /* Gauss elimination */
     
     for (k=0; k < size; k++) {
         if ( fabs(A[k][k]) <= 1e-16 ) {
             printf("Pivot index %d  ",k);
             printf("Pivot value %e  ",A[k][k]);
             Error("Cannot eliminate with such a pivot"); }
         for (i = k+1 ; i <  size; i++) {
             factor = A[i][k] / A[k][k];
             for (j = k+1 ; j < size; j++) 
                 A[i][j] = A[i][j] - A[k][j] * factor;
             B[i] = B[i] - B[k] * factor; }}
     
     /* Back-substitution */
     
     for (i = size-1; i >= 0 ; i--) {
         factor = 0;
         for (j = i+1 ; j < size; j++)
             factor += A[i][j] * B[j];
         B[i] = ( B[i] - factor)/A[i][i]; }
     
     return(mySystem->B);    
 }
 
 void  femFullSystemConstrain(femFullSystem *mySystem, 
                              int myNode, double myValue) 
 {
     double  **A, *B;
     int     i, size;
     
     A    = mySystem->A;
     B    = mySystem->B;
     size = mySystem->size;
     
     for (i=0; i < size; i++) {
         B[i] -= myValue * A[i][myNode];
         A[i][myNode] = 0; }
     
     for (i=0; i < size; i++) 
         A[myNode][i] = 0; 
     
     A[myNode][myNode] = 1;
     B[myNode] = myValue;
 }
 
 
 femProblem *femElasticityCreate(femGeo* theGeometry, 
                   double E, double nu, double rho, double g, femElasticCase iCase)
 {
     femProblem *theProblem = malloc(sizeof(femProblem));
     theProblem->E   = E;
     theProblem->nu  = nu;
     theProblem->g   = g;
     theProblem->rho = rho;
     
     if (iCase == PLANAR_STRESS) {
         theProblem->A = E/(1-nu*nu);
         theProblem->B = E*nu/(1-nu*nu);
         theProblem->C = E/(2*(1+nu)); }
     else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
         theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
         theProblem->B = E*nu/((1+nu)*(1-2*nu));
         theProblem->C = E/(2*(1+nu)); }
 
     theProblem->planarStrainStress = iCase;
     theProblem->nBoundaryConditions = 0;
     theProblem->conditions = NULL;
     
     int size = 2*theGeometry->theNodes->nNodes;
     theProblem->constrainedNodes = malloc(size*sizeof(int));
     theProblem->soluce = malloc(size*sizeof(double));
     theProblem->residuals = malloc(size*sizeof(double));
     for (int i=0; i < size; i++) {
         theProblem->constrainedNodes[i] = -1;
         theProblem->soluce[i] = 0.0;
         theProblem->residuals[i] = 0.0;}
 
 
     
     theProblem->geometry = theGeometry;  
     if (theGeometry->theElements->nLocalNode == 3) {
         theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
         theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
     if (theGeometry->theElements->nLocalNode == 4) {
         theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
         theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
     theProblem->spaceEdge    = femDiscreteCreate(2,FEM_EDGE);
     theProblem->ruleEdge     = femIntegrationCreate(2,FEM_EDGE); 
     theProblem->system       = femFullSystemCreate(size); 
 
     femDiscretePrint(theProblem->space);   
     femDiscretePrint(theProblem->spaceEdge);  
   
     return theProblem;
 }
 
 void femElasticityFree(femProblem *theProblem)
 {
     femFullSystemFree(theProblem->system);
     femIntegrationFree(theProblem->rule);
     femDiscreteFree(theProblem->space);
     femIntegrationFree(theProblem->ruleEdge);
     femDiscreteFree(theProblem->spaceEdge);
     free(theProblem->conditions);
     free(theProblem->constrainedNodes);
     free(theProblem->soluce);
     free(theProblem->residuals);
     free(theProblem);
 }
     

 void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value)
 {
     int iDomain = geoGetDomain(nameDomain);
     // Suppression des vérifs NULL pour concision (elles sont dans Error() ou implicites)
     if (iDomain == -1)  Error("Undefined domain :-(");

     femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
     // Ajout vérif malloc (important)
     if (!theBoundary) Error("AddBoundaryCondition: malloc failed for femBoundaryCondition.");
     theBoundary->domain = theProblem->geometry->theDomains[iDomain];
     theBoundary->value = value;
     theBoundary->type = type;
     theProblem->nBoundaryConditions++;
     int size = theProblem->nBoundaryConditions;
     int indexNewCondition = size - 1; // Index de la condition actuelle

     // Utilisation de realloc pour le tableau de conditions
     femBoundaryCondition **newConditions = realloc(theProblem->conditions, size * sizeof(femBoundaryCondition*));
      if (!newConditions) {
          theProblem->nBoundaryConditions--; free(theBoundary); // Annule et libère
          Error("AddBoundaryCondition: realloc failed for conditions array."); return; }
     theProblem->conditions = newConditions;
     theProblem->conditions[indexNewCondition] = theBoundary;

     // Appliquer contraintes seulement pour Dirichlet
     int shift=-1;
     if (type == DIRICHLET_X)  shift = 0;
     if (type == DIRICHLET_Y)  shift = 1;
     if (shift == -1) return; // Fin si ce n'est pas Dirichlet

     // Préparation pour parcourir les noeuds du domaine
     femDomain *theDomain = theBoundary->domain;
     femMesh   *theEdgesMesh = theDomain->mesh;
     int *domainEdges = theDomain->elem;
     int nEdgesInDomain = theDomain->nElem;
     int nTotalNodes = theProblem->geometry->theNodes->nNodes; // Nb noeuds APRES fix

     // Boucle sur les arêtes du domaine
     for (int e=0; e<nEdgesInDomain; e++) {
         int edgeIndex = domainEdges[e];
          // Vérification minimale de l'indice d'arête
          if (edgeIndex < 0 || edgeIndex >= theEdgesMesh->nElem) continue;

         // Boucle sur les 2 noeuds de l'arête
         for (int i=0; i<theEdgesMesh->nLocalNode; i++) { // nLocalNode == 2
             int nodeIndex = theEdgesMesh->elem[edgeIndex * theEdgesMesh->nLocalNode + i];

             // ================== CHANGEMENT NÉCESSAIRE ==================
             // Ignorer le nœud s'il a été supprimé par geoMeshFix
             if (nodeIndex == -1) {
                 continue; // Passe au nœud suivant
             }
             // ============================================================

             // Vérification que l'indice (non -1) est dans les bornes valides
             if (nodeIndex >= 0 && nodeIndex < nTotalNodes) {
                 int dofIndex = 2 * nodeIndex + shift;
                  // Optionnel : Avertissement si écrasement (conservé car utile)
                 if (theProblem->constrainedNodes[dofIndex] != -1 && theProblem->constrainedNodes[dofIndex] != indexNewCondition) {
                    fprintf(stderr,"Warning: AddBC: Node %d (DOF %d) already constrained. Overwriting.\n", nodeIndex, dofIndex);
                 }
                 // Appliquer la contrainte
                 theProblem->constrainedNodes[dofIndex] = indexNewCondition;
             } else {
                  // Log l'erreur si l'indice est invalide (autre que -1)
                  fprintf(stderr, "Warning: AddBC: Invalid node index %d (max: %d) found for edge %d.\n", nodeIndex, nTotalNodes-1, edgeIndex);
             }
         } // Fin boucle noeuds
     } // Fin boucle arêtes
 }



 
 void femElasticityPrint(femProblem *theProblem)  
 {    
     printf("\n\n ======================================================================================= \n\n");
     printf(" Linear elasticity problem \n");
     printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
     printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
     printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
     printf("   Gravity         g   = %14.7e [m/s2]\n",theProblem->g);
     
     if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
     if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
     if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");
 
     printf("   Boundary conditions : \n");
     for(int i=0; i < theProblem->nBoundaryConditions; i++) {
           femBoundaryCondition *theCondition = theProblem->conditions[i];
           double value = theCondition->value;
           printf("  %20s :",theCondition->domain->name);
           if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
           if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value); 
           if (theCondition->type==NEUMANN_X)    printf(" imposing %9.2e as the horizontal force desnity \n",value); 
           if (theCondition->type==NEUMANN_Y)    printf(" imposing %9.2e as the vertical force density \n",value);}
     printf(" ======================================================================================= \n\n");
 }
 
 double femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y)){
     femIntegration *theRule = theProblem->rule;
     femGeo         *theGeometry = theProblem->geometry;
     femNodes       *theNodes = theGeometry->theNodes;
     femMesh        *theMesh = theGeometry->theElements;
     femDiscrete    *theSpace = theProblem->space;
 
     double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
     int iElem,iInteg,i,map[4];
     int nLocal = theMesh->nLocalNode;
     double value = 0.0;
     for (iElem = 0; iElem < theMesh->nElem; iElem++) {
         for (i=0; i < nLocal; i++) {
             map[i]  = theMesh->elem[iElem*nLocal+i];
             x[i]    = theNodes->X[map[i]];
             y[i]    = theNodes->Y[map[i]];} 
         for (iInteg=0; iInteg < theRule->n; iInteg++) {    
             double xsi    = theRule->xsi[iInteg];
             double eta    = theRule->eta[iInteg];
             double weight = theRule->weight[iInteg];  
             femDiscretePhi2(theProblem->space,xsi,eta,phi);
             femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
             double dxdxsi = 0.0;
             double dxdeta = 0.0;
             double dydxsi = 0.0; 
             double dydeta = 0.0;
             for (i = 0; i < theSpace->n; i++) {  
                 dxdxsi += x[i]*dphidxsi[i];       
                 dxdeta += x[i]*dphideta[i];   
                 dydxsi += y[i]*dphidxsi[i];   
                 dydeta += y[i]*dphideta[i]; }
             double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
             for (i = 0; i < theProblem->space->n; i++) {    
                 value += phi[i] * f(x[i],y[i]) * jac * weight; }}}
     return value;
 
 } 
 
 
 
 double femMin(double *x, int n) 
 {
     double myMin = x[0];
     int i;
     for (i=1 ;i < n; i++) 
         myMin = fmin(myMin,x[i]);
     return myMin;
 }
 
 double femMax(double *x, int n) 
 {
     double myMax = x[0];
     int i;
     for (i=1 ;i < n; i++) 
         myMax = fmax(myMax,x[i]);
     return myMax;
 }
 
 
 void femError(char *text, int line, char *file)                                  
 { 
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
     exit(69);                                                 
 }
 
 void femErrorGmsh(int ierr, int line, char *file)                                  
 { 
     if (ierr == 0)  return;
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
     gmshFinalize(NULL);                                        
     exit(69);                                                 
 }
 
 void femErrorScan(int test, int line, char *file)                                  
 { 
     if (test >= 0)  return;
     
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
     exit(69);                                       
 }
 
 void femWarning(char *text, int line, char *file)                                  
 { 
     printf("\n-------------------------------------------------------------------------------- ");
     printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
     printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
 }


 double femSolverGet(femSolver *mySolver,int i,int j)
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        //case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        //case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]); 
}