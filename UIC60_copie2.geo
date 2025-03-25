SetFactory("OpenCASCADE");

// ========================
// \+ Définition des Points
// ========================
Point(1) = {-75.0, 0, 0, 1.0};              // Extrémité gauche semelle
Point(2) = {75.0, 0, 0, 1.0};               // Extrémité droite semelle
Point(3) = {75.0, 11.5, 0, 1.0};            // Début inclinaison semelle droite
Point(4) = {50.0, 15.5, 0, 1.0};            // Fin inclinaison semelle droite
Point(5) = {12.37, 27.5, 0, 1.0};           // Nouveau point de transition
Point(6) = {-12.37, 27.5, 0, 1.0};          // Symétrique à gauche
Point(7) = {-50.0, 15.5, 0, 1.0};           // Fin inclinaison semelle gauche
Point(8) = {-75.0, 11.5, 0, 1.0};           // Début inclinaison semelle gauche
Point(10) = {8.25, 31.5, 0, 1.0};           // Début âme droite
Point(11) = {8.25, 121, 0, 1.0};            // Âme droite (début transition champignon)
Point(12) = {-8.25, 121, 0, 1.0};           // Âme gauche (début transition champignon)
Point(13) = {-8.25, 31.5, 0, 1.0};          // Début âme gauche
Point(14) = {37.85, 134.3, 0, 1.0};         // Extrémité droite champignon
Point(15) = {37.85, 157.7, 0, 1.0};         // Point intermédiaire droit
Point(16) = {0.0, 172.0, 0, 1.0};           // Centre du rail (sommet champignon)
Point(17) = {-37.85, 157.7, 0, 1.0};        // Point intermédiaire gauche
Point(18) = {-37.85, 134.3, 0, 1.0};        // Extrémité gauche champignon
Point(19) = {12.37, 125, 0, 1.0};           // Âme droite (fin transition champignon)
Point(20) = {-12.37, 125, 0, 1.0};          // Âme gauche (fin transition champignon)
Point(21) = {26.7, 170, 0, 1.0};            // Extrémité droite2 champignon
Point(22) = {-26.7, 170, 0, 1.0};           // Extrémité gauche2 champignon
Point(23) = {10.5, 172, 0, 1.0};            // Extrémité droite3 champignon
Point(24) = {-10.5, 172, 0, 1.0};           // Extrémité gauche3 champignon
Point(25) = {8.9, 29.1, 0, 1.0};
Point(26) = {-33.8, 164.3, 0, 1.0};
Point(27) = {-19.4, 171, 0, 1.0};
Point(28) = {33.8, 164.3, 0, 1.0};
Point(29) = {19.4, 171, 0, 1.0};
Point(30) = {-8.9, 29.1, 0, 1.0};

// ===========================================
// \+ Définition des Courbes (Lines & Splines)
// ===========================================
Line(1) = {1, 2};               // Base semelle
Line(2) = {2, 3};               // Inclinaison semelle droite
Line(3) = {3, 4};               // Transition semelle droite
Line(4) = {4, 5};               // Connexion à l'âme droite
Line(5) = {6, 7};               // Connexion à l'âme gauche
Line(6) = {7, 8};               // Inclinaison semelle gauche
Line(7) = {8, 1};               // Fermeture boucle semelle gauche
Line(9) = {10, 11};             // Âme droite
Line(10) = {12, 13};            // Âme gauche
Line(11) = {11, 19};            // liaison âme droite
Line(12) = {20, 12};            // liaison âme gauche
Line(13) = {19, 14};            // Arrondi champignon droit
Line(14) = {14, 15};            // Partie intermédiaire droite
Line(15) = {17, 18};            // Partie intermédiaire gauche
Line(16) = {18, 20};            // Partie intermédiaire gauche
Line(17) = {23, 16};            // Arrondi champignon droite
Line(18) = {16, 24};            // Arrondi champignon gauche
Spline(34) = {5, 25, 10};
Spline(35) = {22, 26, 17};
Spline(36) = {24, 27, 22};
Spline(37) = {21, 29, 23};
Spline(38) = {15, 28, 21};
Spline(39) = {13, 30, 6};

// ========================
// \+ 📌 SURFACES
// ========================
Curve Loop(1) = {1, 2, 3, 4, 34, 9, 11, 13, 14, 38, 37, 17, 18, 36, 35, 15, 16, 12, 10, 39, 5, 6, 7};  
Plane Surface(1) = {1};

// ========================
// \+ 📌 PHYSICAL GROUPS
// ========================
//Physical Surface("Rail") = {1};                                     // Référence à Plane Surface(1)
//Physical Curve("Bottom") = {39, 5, 6, 7, 1, 2, 3, 4, 34};           // Référence à des courbes
//Physical Curve("Patine") = {17, 18};                                // Référence à des courbes
//Physical Curve("Gauche") = {12, 10};                                // Référence à des courbes
//Physical Curve("Droite") = {11, 9};                                 // Référence à des courbes