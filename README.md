# LEPL1110 | FEM Rail UIC60
---

## Prérequis

Pour compiler et exécuter ce projet, vous aurez besoin de :

1.  **Un compilateur C** (ex: GCC, Clang)
2.  **CMake** (version 3.10 ou supérieure recommandée)
3.  **GMSH** (la bibliothèque `fem.c` l'utilise pour générer le maillage à partir de `UIC60.geo`)
4.  **Python 3** (pour le script d'analyse `load_plot.py`)
5.  **Bibliothèques Python:** `matplotlib` et `numpy`. Installez-les via pip :
    ```bash
    pip install matplotlib numpy
    ```
6.  **Bibliothèques OpenGL/GLFW:** Nécessaires pour la partie visualisation `glfem`. CMake devrait tenter de les trouver. L'installation dépend de votre système (ex: `sudo apt-get install libglfw3-dev libgles2-mesa-dev` sur Debian/Ubuntu).

---

## Compilation

Le projet utilise CMake pour la configuration de la compilation. Suivez ces étapes depuis le répertoire racine du projet (`RAILFEM/`) :

1.  **Créer un répertoire de build :**
    ```bash
    mkdir build
    ```

2.  **Naviguer dans le répertoire de build :**
    ```bash
    cd build
    ```

3.  **Configurer le projet avec CMake :**
    ```bash
    cmake ..
    ```
    *(Assurez-vous que CMake trouve toutes les dépendances nécessaires, notamment GMSH et GLFW/OpenGL si la visualisation est activée).*

4.  **Compiler le projet :**
    ```bash
    make
    ```
    Cela devrait générer un exécutable nommé `myFem` (ou un autre nom si défini différemment dans `CMakeLists.txt`) dans le répertoire `build/`.

---

## Exécution de la Simulation

Une fois la compilation réussie :

1.  **Restez dans le répertoire `build/`** (ou naviguez-y : `cd build`).
2.  **Lancez l'exécutable :**
    ```bash
    ./myFem
    ```

L'exécution du programme va :
*   Appeler GMSH pour générer le maillage basé sur `Projet-EF/UIC60.geo`.
*   Importer et potentiellement corriger le maillage (`geoMeshFix`).
*   Définir le problème d'élasticité (matériau, CL).
*   Résoudre le système par éléments finis.
*   Sauvegarder le maillage final dans `data/mesh.txt`.
*   Sauvegarder les résultats numériques (déplacements, etc.) dans `data/elasticity.txt`.
*   Ouvrir une fenêtre OpenGL pour la visualisation interactive des résultats (si activée et fonctionnelle). Utilisez les touches indiquées dans la console (V, D, X, Y, N, ESC) pour interagir.

---

## Génération des Graphiques d'Analyse

Le script Python `load_plot.py` génère les graphiques montrant la relation charge-déplacement au point spécifique (19.4, 171.0). Les données pour ces graphiques sont actuellement *codées en dur* dans le script.

1.  **Depuis le répertoire racine du projet (`RAILFEM/`)**, exécutez le script :
    ```bash
    python Projet-EF/analyse/load_plot.py
    ```
    *(Alternativement, naviguez dans `Projet-EF/analyse/` et exécutez `python load_plot.py`)*

2.  Cela va :
    *   Afficher les deux graphiques à l'écran.
    *   Sauvegarder les graphiques sous forme de fichiers PNG dans le répertoire où le script est exécuté (par défaut, dans `Projet-EF/analyse/`):
        *   `rail_load_vs_displacement_P29.png` (toutes les charges)
        *   `rail_load_vs_displacement_P29_limited.png` (charges <= 1e8 N/m)

---

## Fichiers de Sortie Principaux

*   `build/myFem`: Exécutable principal de la simulation.
*   `data/mesh.txt`: Maillage utilisé pour la simulation (nœuds, éléments).
*   `data/elasticity.txt`: Résultats numériques (peut contenir les déplacements nodaux).
*   `Projet-EF/analyse/rail_load_vs_displacement_P29.png`: Graphique charge-déplacement (toutes charges).
*   `Projet-EF/analyse/rail_load_vs_displacement_P29_limited.png`: Graphique charge-déplacement (charges limitées).

---

## Auteurs

*   Théodore Moulaert
*   Edouard Meurant