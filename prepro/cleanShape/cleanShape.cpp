// Copyright (C) shapeSurvey <vincent.richefeu@3sr-grenoble.fr>
//
// This file is part of mbox.
//
// shapeSurvey can not be copied and/or distributed without the express
// permission of the authors.
// It is coded for academic purposes.
//
// Note
// Without a license, the code is copyrighted by default.
// People can read the code, but they have no legal right to use it.
// To use the code, you must contact the author directly and ask permission.
//
// Convert STL to shapeSurvey format
// Author: Vincent.Richefeu@3sr-grenoble.fr
// Lab 3SR, Grenoble
// 2016

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>

#include <map>
#include <set>
#include <vector>

#include <sstream>
#include <stdint.h>
#include <string>

#define ZERO 1e-12

#include "AABB.hpp"
#include "CmdLine.h"
#include "Shape.hpp"
#include "fileTool.hpp"
#include "message.hpp"
#include "vec3.hpp"

std::vector<Shape> Shapes;

auto isTowFacesConnected = [](const std::vector<size_t> &face1, const std::vector<size_t> &face2) {
    size_t nb_common_vertex = 0;
    for (size_t i = 0; i < face1.size(); i++) {
        for (size_t j = 0; j < face2.size(); j++) {
            if (face1[i] == face2[j])
                nb_common_vertex++;
        }
    }
    return (nb_common_vertex >= 2);
};

auto commonEdgeFacesConnected = [](const std::vector<size_t> &face1, const std::vector<size_t> &face2) {
    std::pair<size_t, size_t> common_edge;
    size_t nb_common_vertex = 0;

    for (size_t i = 0; i < face1.size(); i++) {
        for (size_t j = 0; j < face2.size(); j++) {
            if (face1[i] == face2[j]) {
                nb_common_vertex++;
                if (nb_common_vertex == 1)
                    common_edge.first = face2[j];
                if (nb_common_vertex == 2)
                    common_edge.second = face2[j];
            }
        }
    }
    return common_edge;
};

void removeEgde(std::vector<std::pair<size_t, size_t>> &vedge, std::pair<size_t, size_t> edge) {
    std::vector<std::pair<size_t, size_t>>::iterator it = vedge.begin();
    while (it != vedge.end()) {
        std::pair<size_t, size_t> edg = (std::pair<size_t, size_t>)*it;
        if (((edg.first == edge.first) && (edg.second == edge.second)) || ((edg.first == edge.second) && (edg.second == edge.first))) {
            vedge.erase(it);
            break;
        } else
            ++it;
    }
}

void removeFace(std::vector<std::vector<size_t>> &face, size_t f1) {
    // erase the f1 th element
    face.erase(face.begin() + f1);
}

void addEgde(std::vector<std::pair<size_t, size_t>> &vedge, std::pair<size_t, size_t> edge) {
    for (size_t i = 0; i < vedge.size(); i++) {
        std::pair<size_t, size_t> edg = vedge[i];
        if (((edg.first == edge.first) && (edg.second == edge.second)) || ((edg.first == edge.second) && (edg.second == edge.first)))
            return;
    }
    vedge.push_back(edge);
}

void concatener_deux_faces(std::vector<std::vector<size_t>> &face, size_t fmax, size_t fmin, std::pair<size_t, size_t> edge) {
    //edge contient l'orientation de face max

    int debut = 0;
    int fin = 0;
    //il me faut trouver la position de edge.first dans fmin
    for (size_t i = 0; i < face[fmin].size(); i++) {
        if (edge.first == face[fmin][i])
            debut = i;
        if (edge.second == face[fmin][i])
            fin = i;
    }

    std::vector<int> append;

    int nbvmax = face[fmin].size() - 2; //nombre de vertex maximal a ajouter

    append.resize(nbvmax);

    if (fin - debut > 0) //le cas edial les deux se suivent
    {
        if (fin - debut == 1) {
            int nbv = 0;
            for (int i = debut - 1; i >= 0; i--) {
                append[nbv] = face[fmin][i];
                nbv++;
            }

            for (int i = face[fmin].size() - 1; i > fin; i--) {
                append[nbv] = face[fmin][i];
                nbv++;
            }
        } else //extrimite
        {
            int nbv = 0;
            for (int i = debut + 1; i < fin; i++) {
                append[nbv] = face[fmin][i];
                nbv++;
            }
        }
    } else //le cas les deux se suivnet mais il faut inverser
    {
        if (fin - debut == -1) {
            int nbv = 0;
            for (size_t i = debut + 1; i < face[fmin].size(); i++) {
                append[nbv] = face[fmin][i];
                nbv++;
            }

            for (int i = 0; i < fin; i++) {
                append[nbv] = face[fmin][i];
                nbv++;
            }
        } else //extrimite
        {
            int nbv = 0;
            for (int i = debut - 1; i > fin; i--) {
                append[nbv] = face[fmin][i];
                nbv++;
            }
        }
    }

    std::vector<int> ids;

    //il me faut trouver la position de edge.first dans fmin
    for (size_t i = 0; i < face[fmax].size(); i++) {
        if (edge.first == face[fmax][i])
            debut = i;
        if (edge.second == face[fmax][i])
            fin = i;
    }

    if (fin - debut == 1) {
        //les deux se suivent
        for (int i = 0; i <= debut; i++) {
            ids.push_back(face[fmax][i]);
        }
        for (size_t i = 0; i < append.size(); i++) {
            ids.push_back(append[i]);
        }
        for (size_t i = fin; i < face[fmax].size(); i++) {
            ids.push_back(face[fmax][i]);
        }
    } else {
        //les deux sont à l'extremmité
        for (size_t i = 0; i < face[fmax].size(); i++) {
            ids.push_back(face[fmax][i]);
        }
        for (size_t i = 0; i < append.size(); i++) {
            ids.push_back(append[i]);
        }
    }

    face[fmax].clear();
    for (size_t i = 0; i < ids.size(); i++) {
        face[fmax].push_back(ids[i]);
    }

    //j enleve la face
    removeFace(face, fmin);
}

vec3r get_normal_face(size_t &f, Shape &shape) {
    vec3r n;
    std::vector<size_t> face = shape.face[f];

    //je cherche les deux vecteur qui ont le plus grand produit scalaire :
    std::vector<vec3r> vecteurs;

    for (size_t i = 0; i < face.size(); i++) {
        vec3r vi = shape.vertex[face[i]];
        for (size_t j = i + 1; j < face.size(); j++) {
            vec3r vj = shape.vertex[face[j]];
            vec3r vji = vj - vi;
            vecteurs.push_back(vji);
        }
    }

    double maxdot = 0.;
    size_t maxi = 0, maxj = 0;

    for (size_t i = 0; i < vecteurs.size(); i++) {
        vec3r vij_i = vecteurs[i];
        for (size_t j = i + 1; j < vecteurs.size(); j++) {
            vec3r vij_j = vecteurs[j];

            if (fabs(dot(vij_i, vij_j)) > maxdot) {
                maxdot = fabs(dot(vij_i, vij_j));
                maxi = i;
                maxj = j;
            }
        }
    }
    n = vecteurs[maxi] ^ vecteurs[maxj];
    n.normalized();
    return n;
}

int cleanShape(Shape &shape) {
    double dotscalar = 0.; //au plan

    for (size_t i = 0; i < shape.face.size(); i++) {
        vec3r ni = get_normal_face(i, shape);
        for (size_t j = i + 1; j < shape.face.size(); j++) {
            if (isTowFacesConnected(shape.face[i], shape.face[j])) {
                vec3r nj = get_normal_face(j, shape);
                dotscalar = fabs(fabs(dot(ni, nj)) - 1.);
                if (dotscalar <= ZERO) {
                    //là il faut penser à une methode de fusion des faces et enlever l'edge commun

                    //decider quel face il faut suprimer d abord
                    //ici c'est la face qui est en deuxieme parametre
                    std::pair<size_t, size_t> edge;

                    if (shape.face[i].size() > shape.face[j].size())
                        edge = commonEdgeFacesConnected(shape.face[i], shape.face[j]);
                    else
                        edge = commonEdgeFacesConnected(shape.face[j], shape.face[i]);

                    removeEgde(shape.edge, edge);

                    //concatener les deux faces

                    if (shape.face[i].size() > shape.face[j].size()) {
                        concatener_deux_faces(shape.face, i, j, edge);
                    } else {
                        concatener_deux_faces(shape.face, j, i, edge);
                    }
                    return 1;
                }
            }
        }
    }
    return 0;
}

void loadShapes(const char *fileName) {
    // If a library file is in the running folder, so it is preferably used
    std::string ModFileName(fileName);
    std::string LocalFileName = fileTool::GetFileName(ModFileName) + "." + fileTool::GetFileExt(ModFileName);
    if (fileTool::fileExists(LocalFileName.c_str())) {
        ModFileName = LocalFileName;
    }

    if (!fileTool::fileExists(ModFileName.c_str())) {
        std::cerr << msg::warn() << "@CleanShape::loadShapes, shape library named '" << ModFileName << "' has not been found."
                  << msg::normal() << std::endl
                  << std::endl;
        return;
    }

    Shapes.clear();

    std::vector<Shape> shapes;

    std::ifstream is(ModFileName.c_str());

    std::string token;
    is >> token;
    while (is) {
        if (token == "<") {
            Shape S;
            S.read(is);
            shapes.push_back(S);
        }
        is >> token;
    }

    Shapes.resize(shapes.size());
    for (size_t i = 0; i < shapes.size(); i++) {
        Shapes[i] = shapes[i];
    }

    std::cout << "Number of shapes found in the library file " << ModFileName << ": " << shapes.size() << std::endl;
}

int main(int argc, char const *argv[]) {
    std::string shpFileName;

    try {
        TCLAP::CmdLine cmd("Clean a shape file that will be used by Rockable", ' ', "0.3");
        TCLAP::UnlabeledValueArg<std::string> nameArg("input", "Name of the input shape file", true, "file.shp",
                                                      "shp file");
        cmd.add(nameArg);

        cmd.parse(argc, argv);

        shpFileName = nameArg.getValue();
    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
    }

    loadShapes(shpFileName.c_str());

    for (size_t s = 0; s < Shapes.size(); s++) {
        Shape &S = Shapes[s];
        //A utiliser uniquement lorsqu on fais es polyhedres reguliers
        int netoyage = 0;
        do {
            netoyage = cleanShape(S);
        } while (netoyage != 0);
        S.defineVertexConnectivity();
        S.buildOBBtree();
        if (S.preCompDone == 'n') {
            S.massProperties();
        }
    }

    shpFileName = "m_" + shpFileName;

    std::ofstream os(shpFileName.c_str());

    for (size_t s = 0; s < Shapes.size(); s++) {
        Shape &S = Shapes[s];
        S.write(os);
    }
    return 0;
}
