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

Shape pale;
Shape base_cylindre;

auto searchvertex = [](std::vector<vec3r> &arr, vec3r &v) {
    auto it = std::find(arr.begin(), arr.end(), v);
    if (it != arr.end())
        return (int)(it - arr.begin());
    return (int)-1;
};

auto searchedge = [](std::vector<std::pair<size_t, size_t>> &arr, std::pair<size_t, size_t> &e) {
    auto it = std::find(arr.begin(), arr.end(), e);
    if (it != arr.end())
        return (int)(it - arr.begin());
    return (int)-1;
};

auto create = [](vec3r &v0, vec3r &v1, const int &res, const int &max_nb, Shape &shape) {
    std::vector<vec3r> &vertexes = shape.vertex;
    std::vector<std::pair<size_t, size_t>> &edges = shape.edge;
    std::vector<std::vector<size_t>> &faces = shape.face;

    double dtheta = 2. * M_PI / (double)res;
    for (int i = 0; i < max_nb; i++) {

        vec3r v[4];

        v[0].x = v0.x * cos(dtheta * ((i) % res));
        v[0].y = v0.x * sin(dtheta * ((i) % res));
        v[0].z = v0.z;

        v[1].x = v1.x * cos(dtheta * ((i) % res));
        v[1].y = v1.x * sin(dtheta * ((i) % res));
        v[1].z = v1.z;

        v[2].x = v1.x * cos(dtheta * ((i + 1) % res));
        v[2].y = v1.x * sin(dtheta * ((i + 1) % res));
        v[2].z = v1.z;

        v[3].x = v0.x * cos(dtheta * ((i + 1) % res));
        v[3].y = v0.x * sin(dtheta * ((i + 1) % res));
        v[3].z = v0.z;

        size_t iadd[4];

        //vertex;
        for (size_t vi = 0; vi < 4; vi++) {

            int num = searchvertex(vertexes, v[vi]);

            if (num == -1) {
                iadd[vi] = vertexes.size();
                vertexes.push_back(v[vi]);
            } else {
                iadd[vi] = (size_t)num;
            }
        }
        //edges
        std::pair<size_t, size_t> e[4];

        e[0].first = iadd[0];
        e[0].second = iadd[1];

        e[1].first = iadd[1];
        e[1].second = iadd[2];

        e[2].first = iadd[2];
        e[2].second = iadd[3];

        e[3].first = iadd[3];
        e[3].second = iadd[0];

        for (size_t ei = 0; ei < 4; ei++) {
            if (e[ei].first < e[ei].second)
                std::swap(e[ei].first, e[ei].second);
            int num = searchedge(edges, e[ei]);
            if (num == -1)
                edges.push_back(e[ei]);
        }

        //faces
        std::vector<size_t> idf;
        for (size_t vi = 0; vi < 4; vi++)
            idf.push_back(iadd[vi]);
        faces.push_back(idf);
    }
};

void make_arm(Shape &S, const int &resolution) {

    S.name = "bras";

    //construction du bas de la pale
    vec3r p0 = vec3r(1.5, 0., -8.45);
    vec3r p1 = vec3r(2.875, 0., -7.2);
    create(p0, p1, resolution, resolution, S);

    std::vector<size_t> idf;
    for (size_t v = 0; v < S.vertex.size(); v++) {
        if (S.vertex[v].z == p0.z)
            idf.push_back(v);
    }
    S.face.push_back(idf);

    //
    p0 = vec3r(2.875, 0., -7.2);
    p1 = vec3r(2.875, 0., -6.2);
    create(p0, p1, resolution, resolution, S);

    //
    p0 = vec3r(2.875, 0., -6.2);
    p1 = vec3r(5., 0., -6.2);
    create(p0, p1, resolution, resolution, S);

    //
    p0 = vec3r(5., 0., -6.2);
    p1 = vec3r(5., 0., 0.);
    create(p0, p1, resolution, resolution, S);

    //
    p0 = vec3r(5., 0., 0.);
    p1 = vec3r(3.5, 0., 1.5);
    create(p0, p1, resolution, resolution, S);

    //
    p0 = vec3r(3.5, 0., 1.5);
    p1 = vec3r(3.5, 0., 133.5);
    create(p0, p1, resolution, resolution, S);

    idf.clear();
    for (size_t v = 0; v < S.vertex.size(); v++) {
        if (S.vertex[v].z == p1.z)
            idf.push_back(v);
    }

    S.face.push_back(idf);

    S.defineVertexConnectivity();
    S.buildOBBtree();

    std::string shpbras = "bras.shp";
    std::ofstream os_bras(shpbras.c_str());
    S.write(os_bras);
}

void makeBaseCylindre(Shape &S, const int &resolution) {

    S.name = "baseCyl";

    //construction du bas de la pale
    vec3r p0 = vec3r(25, 0., -17.);
    vec3r p1 = vec3r(14, 0., -15.);
    create(p0, p1, resolution, resolution, S);

    std::vector<size_t> idf;
    for (size_t v = 0; v < S.vertex.size(); v++) {
        if (S.vertex[v].z == p0.z)
            idf.push_back(v);
    }
    S.face.push_back(idf);

    idf.clear();
    for (size_t v = 0; v < S.vertex.size(); v++) {
        if (S.vertex[v].z == p1.z)
            idf.push_back(v);
    }

    S.face.push_back(idf);

    S.defineVertexConnectivity();
    S.buildOBBtree();
}

void make_pale(Shape &S, std::string &fileName) {

    size_t nv = S.vertex.size();
    std::vector<vec3r> vertex = S.vertex;
    std::vector<std::pair<size_t, size_t>> edges = S.edge;
    std::vector<std::vector<size_t>> faces = S.face;

    S = Shape();

    std::ifstream is(fileName.c_str());
    std::string token;
    is >> token;
    while (is) {
        if (token == "<") {
            S.read(is);
        }
        is >> token;
    }

    //vertex
    for (size_t i = 0; i < S.vertex.size(); i++) {
        vertex.push_back(S.vertex[i]);
    }
    //edges
    for (size_t i = 0; i < S.edge.size(); i++) {
        std::pair<size_t, size_t> e = S.edge[i];
        e.first += nv;
        e.second += nv;
        edges.push_back(e);
    }
    //faces
    for (size_t i = 0; i < S.face.size(); i++) {
        std::vector<size_t> ff = S.face[i];

        for (size_t f = 0; f < ff.size(); f++)
            ff[f] += nv;
        faces.push_back(ff);
    }

    S.vertex = vertex;
    S.edge = edges;
    S.face = faces;

    S.defineVertexConnectivity();
    S.buildOBBtree();
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

    int resolution = 50;

    Shape &S0 = pale;
    pale.radius = 1e-3;
    make_arm(S0, resolution);
    make_pale(S0, shpFileName);
    Shape &S1 = base_cylindre;
    makeBaseCylindre(S1, resolution);

    //write shapes
    shpFileName = "m_" + shpFileName;
    std::ofstream os(shpFileName.c_str());
    S0.write(os);
    S1.write(os);
    os.close();

    return 0;
}