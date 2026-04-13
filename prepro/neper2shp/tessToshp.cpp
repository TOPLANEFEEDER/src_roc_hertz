//
//  tessToshp.cpp
//  rockable
//
//  Created by amarsid on 25/09/2019.
//  Copyright © 2019 amarsid. All rights reserved.
//

#include "tessToshp.hpp"
#include "polyhTool.hpp"
#include <cstdlib>

#define EPSILON 1.e-10

//convert a shape to a poly
int PolyToShape(const poly &block, Shape &shape) {
    shape.vertex.clear();
    shape.edge.clear();
    shape.face.clear();
    size_t nv = block.nodes.size();
    size_t ne = block.edges.size();
    size_t nf = block.faces.size();
    //vetexes
    for (size_t iv = 0; iv < nv; iv++) {
        vec3r v = vec3r();
        v = block.nodes[iv].pos;
        shape.vertex.push_back(v);
    }
    //edges
    for (size_t ie = 0; ie < ne; ie++) {
        std::pair<size_t, size_t> e;
        e.first = block.edges[ie].node0;
        e.second = block.edges[ie].node1;
        shape.edge.push_back(e);
    }
    //faces
    for (size_t ifc = 0; ifc < nf; ifc++) {
        pface f = block.faces[ifc];
        std::vector<size_t> ids;
        for (auto it = f.nodes.begin(); it != f.nodes.end(); it++) {
            size_t node_face = *it;
            ids.push_back(node_face);
        }
        shape.face.push_back(ids);
    }
    return 1;
}

poly cube(vec3r min_v, vec3r max_v) {
    poly p;

    p.nodes.resize(8);
    p.edges.resize(12);
    p.faces.resize(6);

    p.nodes[0].pos = vec3r(min_v.x, min_v.y, min_v.z);
    p.nodes[1].pos = vec3r(max_v.x, min_v.y, min_v.z);
    p.nodes[2].pos = vec3r(max_v.x, min_v.y, max_v.z);
    p.nodes[3].pos = vec3r(min_v.x, min_v.y, max_v.z);
    p.nodes[4].pos = vec3r(min_v.x, max_v.y, min_v.z);
    p.nodes[5].pos = vec3r(max_v.x, max_v.y, min_v.z);
    p.nodes[6].pos = vec3r(max_v.x, max_v.y, max_v.z);
    p.nodes[7].pos = vec3r(min_v.x, max_v.y, max_v.z);

    p.edges[0].node0 = 0;
    p.edges[0].node1 = 1;
    p.edges[1].node0 = 1;
    p.edges[1].node1 = 2;
    p.edges[2].node0 = 2;
    p.edges[2].node1 = 3;
    p.edges[3].node0 = 3;
    p.edges[3].node1 = 0;
    p.edges[4].node0 = 4;
    p.edges[4].node1 = 5;
    p.edges[5].node0 = 5;
    p.edges[5].node1 = 6;
    p.edges[6].node0 = 6;
    p.edges[6].node1 = 7;
    p.edges[7].node0 = 7;
    p.edges[7].node1 = 4;
    p.edges[8].node0 = 0;
    p.edges[8].node1 = 4;
    p.edges[9].node0 = 5;
    p.edges[9].node1 = 1;
    p.edges[10].node0 = 2;
    p.edges[10].node1 = 6;
    p.edges[11].node0 = 3;
    p.edges[11].node1 = 7;

    p.faces[0].nodes = {0, 1, 2, 3};
    p.faces[1].nodes = {4, 5, 6, 7};
    p.faces[2].nodes = {0, 1, 5, 4};
    p.faces[3].nodes = {1, 5, 6, 2};
    p.faces[4].nodes = {2, 6, 7, 3};
    p.faces[5].nodes = {0, 4, 7, 3};

    return p;
}

void tessToshp::addPlan(double theta, double psi) {
    double eps = 1.e-15;
    plan_ plan;
    double d = 0.;
    vec3r n;
    vec3r pos;

    n.x = cos(psi) * sin(theta);
    n.y = sin(psi) * sin(theta);
    n.z = cos(theta);

    pos = R * n;
    d = pos * n;

    plan.pos = pos;
    plan.normal = n;
    plan.d = d;

    bool exist = false;

    for (size_t ip = 0; ip < plans.size(); ip++) {
        plan_ p = plans[ip];
        if (p.d == plan.d) {
            vec3r n = plan.normal - p.normal;
            if (fabs(n.x) < eps && fabs(n.y) < eps && fabs(n.z) < eps) {
                exist = true;
                continue;
            }
        }
    }

    if (!exist)
        plans.push_back(plan);
}

void tessToshp::createShpere() {
    plans.clear();
    double d_theta = M_PI / (double)nbP;
    double d_psi = 2. * M_PI / (double)nbP;

    double theta = 0.;
    double psi = 0.;

    addPlan(theta, psi);

    theta = M_PI;
    psi = 0.;

    addPlan(theta, psi);

    for (int i_psi = 0; i_psi < nbP; i_psi++) {
        for (int i_theta = 1; i_theta < nbP; i_theta++) {
            theta = i_theta * d_theta;
            psi = i_psi * d_psi;
            addPlan(theta, psi);
        }
    }

    std::string plans_str = "./plans";
    std::ofstream out(plans_str.c_str());
    out << std::scientific << std::setprecision(13);
    out << plans.size() << std::endl;
    for (size_t i_plan = 0; i_plan < plans.size(); i_plan++) {
        plan_ p = plans[i_plan];
        out << p.d << '\t' << p.normal.x << '\t' << p.normal.y << '\t' << p.normal.z << std::endl;
    }
}

void tessToshp::createTess() {
    std::string command_tess;

    command_tess = "neper -T -n " + std::to_string(number_of_cells) + " -dim 3 -domain 'planes('./plans')' -morpho centroidal -morphooptialgo 'lloyd(1)' -morphooptistop itermax=500 -o 'out_tess'";
    //command_tess = "neper -T -n " + std::to_string(number_of_cells) + " -dim 3 -o 'out_tess'";
    std::system(command_tess.c_str());
}

void tessToshp::readTess() {
    std::string is = "./out_tess.tess";
    std::ifstream fichier(is.c_str(), std::ios::in);

    if (fichier) {

        std::string line;
        int tmp_int;
        double tmp_double;

        // Cells
        do
            std::getline(fichier, line);
        while (line != std::string(" **cell"));

        fichier >> number_of_cells;
        std::cout << "> number_of_cells " << number_of_cells << std::endl
                  << std::flush;

        // Positions
        do
            std::getline(fichier, line);
        while (line != std::string("  *seed"));
        position_.resize(number_of_cells);
        for (size_t i = 0; i < position_.size(); i++)
            fichier >> position_[i].x >> position_[i].y >> position_[i].z >> tmp_double;

        // Vertices
        do
            std::getline(fichier, line);
        while (line != std::string(" **vertex"));

        fichier >> tmp_int;
        vertex.resize(tmp_int);
        for (size_t i = 0; i < vertex.size(); i++)
            fichier >> vertex[i].id >> vertex[i].v.x >> vertex[i].v.y >> vertex[i].v.z >> tmp_int;

        // Edges
        do
            std::getline(fichier, line);
        while (line != std::string(" **edge"));
        fichier >> tmp_int;
        edge.resize(tmp_int);
        for (size_t i = 0; i < edge.size(); i++)
            fichier >> edge[i].id >> edge[i].ver_1 >> edge[i].ver_2 >> tmp_int;

        // Face
        do
            std::getline(fichier, line);
        while (line != std::string(" **face"));
        fichier >> tmp_int;
        face.resize(tmp_int);
        for (size_t i = 0; i < face.size(); i++) {
            fichier >> face[i].id;
            fichier >> tmp_int;           // number of vertices
            face[i].ver_.resize(tmp_int); //cout << face[i].ver_.size() << std::endl;
            for (size_t j = 0; j < face[i].ver_.size(); j++)
                fichier >> face[i].ver_[j];

            fichier >> tmp_int; // number of edges
            face[i].edge_.resize(tmp_int);
            for (size_t j = 0; j < face[i].edge_.size(); j++)
                fichier >> face[i].edge_[j];

            fichier >> face[i].face_eq_d >> face[i].face_eq_a >> face[i].face_eq_b >> face[i].face_eq_c;
            fichier >> tmp_int >> tmp_int;
            fichier >> tmp_double >> tmp_double >> tmp_double;
        }

        // Polyhedron
        do
            std::getline(fichier, line);
        while (line != std::string(" **polyhedron"));
        fichier >> tmp_int;
        polyhedron.resize(tmp_int);
        for (size_t i = 0; i < polyhedron.size(); i++) {
            fichier >> tmp_int;
            fichier >> tmp_int; // number_of_faces
            polyhedron[i].face_.resize(tmp_int);
            for (size_t j = 0; j < polyhedron[i].face_.size(); j++)
                fichier >> polyhedron[i].face_[j];
        }

        fichier.close();
    } else
        std::cerr << "Impossible d'ouvrir le fichier !" << std::endl;
}

void tessToshp::create_shapes() {
    for (size_t is = 0; is < Shapes.size(); is++) {
        size_t nsub = 0;
        Shape *shp = Shapes[is];
        AABB aabb(shp->vertex);
        poly block = cube(aabb.min, aabb.max);

        poly sub_block;

        for (size_t ip = 0; ip < shape_plans[is].size(); ip++) {
            plan p = shape_plans[is][ip];
            p.normal.normalized();
            p.normal *= -1.;
            p.pos += radius * p.normal;
            //from here, there is intersection
            if (polyhToolCut::cut_poly(block, p, sub_block) == 1) {
                block = sub_block;
                nsub++;
            }
        }
        if (nsub != 0) {
            PolyToShape(block, *shp);
        }
    }

    min_norm_edge = +1.e+20;
    max_norm_edge = -1.e+20;

    for (size_t is = 0; is < Shapes.size(); is++) {
        Shape *shp = Shapes[is];
        for (size_t ie = 0; ie < shp->edge.size(); ie++) {
            vec3r ve = shp->vertex[shp->edge[ie].first] - shp->vertex[shp->edge[ie].second];
            double n_ve = norm(ve);
            min_norm_edge = std::min(min_norm_edge, n_ve);
            max_norm_edge = std::max(max_norm_edge, n_ve);

            if (min_norm_edge == n_ve)
                std::cout << " -Shape number " << is << ", n_ve " << n_ve << std::endl;
        }
    }
}

void tessToshp::extractShapes() {
    size_t np = polyhedron.size();

    Shapes.resize(np);
    face_plans.resize(np);
    shape_plans.resize(np);

    for (size_t ip = 0; ip < polyhedron.size(); ip++) {
        Shape& shp_id = *Shapes[ip];
        vec3r &position = position_[ip];
        std::vector<plan_> &face_plan = face_plans[ip];
        std::vector<plan> &shape_plan = shape_plans[ip];

        shp_id.name = "shape_" + std::to_string(ip);
        shp_id.radius = radius;
        shp_id.volume = 23056;
        shp_id.inertia_mass = vec3r(280.628, 223.074, 128.42);
        shp_id.obb.e[0] = vec3r::unit_x();
        shp_id.obb.e[1] = vec3r::unit_y();
        shp_id.obb.e[2] = vec3r::unit_z();
        shp_id.orientation = quat(0., 0., 0., 1.0);
        shp_id.position = position;
        shp_id.MCnstep = 50000;

        std::vector<vertex_> V;

        // Cherche, classe et supprime les doublons dans la liste des vertex
        for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) { // loop over faces
            int face_id = abs(polyhedron[ip].face_[if1]);
            face_ ff;
            ff.id = face_id;
            int face_index = getIndex(face, ff);
            for (size_t iv1 = 0; iv1 < face[face_index].ver_.size(); iv1++) { // loop over vertex
                vertex_ vv;
                vv.id = face[face_index].ver_[iv1];
                int vertex_index = getIndex(vertex, vv);
                V.push_back(vertex[vertex_index]);
            }
        }

        sort(V.begin(), V.end());                     // Trie le vecteur
        V.erase(unique(V.begin(), V.end()), V.end()); // supprime les doublons

        // Cherche la liste des edges
        std::vector<edge_> E;                                            // Edge id
        for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) { // loop over faces
            int face_id = abs(polyhedron[ip].face_[if1]);
            face_ ff;
            ff.id = face_id;
            int face_index = getIndex(face, ff);
            for (size_t ie1 = 0; ie1 < face[face_index].edge_.size(); ie1++) { // loop over vertex
                edge_ ee;
                ee.id = face[face_index].edge_[ie1];
                int edge_index = getIndex(edge, ee);
                E.push_back(edge[edge_index]);
            }
        }

        sort(E.begin(), E.end());
        E.erase(unique(E.begin(), E.end()), E.end()); // supprime les doublons

        shp_id.vertex.resize(V.size());

        for (size_t iv1 = 0; iv1 < V.size(); iv1++) {
            vec3r v = vec3r();
            v = V[iv1].v;
            shp_id.vertex[iv1] = v;
        }

        AABB aabb(shp_id.vertex);
        vec3r extent = (aabb.max - aabb.min) / 2.;
        vec3r center = (aabb.max + aabb.min) / 2.;
        shp_id.obb.extent = extent;
        shp_id.obb.center = center;

        shp_id.edge.resize(E.size());

        for (size_t ie1 = 0; ie1 < E.size(); ie1++) {
            std::pair<size_t, size_t> e;

            vertex_ vv1;
            vv1.id = E[ie1].ver_1;
            vertex_ vv2;
            vv2.id = E[ie1].ver_2;
            int vertex1_index = getIndex(V, vv1);
            int vertex2_index = getIndex(V, vv2);

            e.first = vertex1_index;
            e.second = vertex2_index;
            shp_id.edge[ie1] = e;
        }

        shp_id.face.resize(polyhedron[ip].face_.size());
        face_plan.resize(shp_id.face.size());
        shape_plan.resize(shp_id.face.size());

        for (size_t if1 = 0; if1 < polyhedron[ip].face_.size(); if1++) {
            // loop over faces
            int face_id = polyhedron[ip].face_[if1];

            face_ ff;
            ff.id = abs(face_id);
            int face_index = getIndex(face, ff);

            shp_id.face[if1].resize(face[face_index].ver_.size());

            face_plan[if1].d = face[face_index].face_eq_d;
            face_plan[if1].normal.x = face[face_index].face_eq_a;
            face_plan[if1].normal.y = face[face_index].face_eq_b;
            face_plan[if1].normal.z = face[face_index].face_eq_c;

            shape_plan[if1].normal = face_plan[if1].normal;
            if (face_id < 0)
                shape_plan[if1].normal *= -1.;

            vec3r pos = vec3r();
            for (size_t iv1 = 0; iv1 < face[face_index].ver_.size(); iv1++) { // loop over vertex
                vertex_ vv;
                vv.id = face[face_index].ver_[iv1];
                int vertex_index = getIndex(V, vv);
                shp_id.face[if1][iv1] = vertex_index;
                pos += V[vertex_index].v;
            }

            int nv = (int)face[face_index].ver_.size();
            if (nv != 0)
                pos /= (double)nv;

            shape_plan[if1].pos = pos;
        }
    }
}

void tessToshp::validate_shape_faces(size_t num_shp) {
    Shape& shp = *Shapes[num_shp];
    std::vector<plan_> shape_plans = face_plans[num_shp];

    //    std::cout<<"______________________________________________"<<std::endl;

    //    std::cout<<"Shape num "<<num_shp<<std::endl;

    double min_edge = +1.e+20;
    double max_edge = -1.e+20;

    for (size_t e = 0; e < shp.edge.size(); e++) {
        vec3r eij = shp.vertex[shp.edge[e].first] - shp.vertex[shp.edge[e].second];
        double d_eij = norm(eij);

        min_edge = std::min(d_eij, min_edge);
        max_edge = std::max(d_eij, max_edge);
    }

    //    std::cout<<"min norm edge "<<min_edge<<std::endl;
    //    std::cout<<"max norm edge "<<max_edge<<std::endl;

    min_norm_edge = std::min(min_norm_edge, min_edge);
    max_norm_edge = std::max(max_norm_edge, max_edge);

    for (size_t f = 0; f < shp.face.size(); f++)
        if (shp.face[f].size() > 3) {
            std::vector<std::pair<size_t, size_t>> edges;

            vec3r n = shape_plans[f].normal;
            double d = shape_plans[f].d;

            //remplir la liste des vertexes de la face
            std::vector<size_t> face_vertexes = shp.face[f];

            std::vector<size_t> out_plan;
            std::vector<size_t> inside_plan;

            double max_res = -1.e+20;
            double min_res = +1.e-20;

            //tester si la face est bien definie
            for (size_t i = 0; i < face_vertexes.size(); i++) {
                double res = n * shp.vertex[face_vertexes[i]] - d;

                max_res = std::max(fabs(res), max_res);
                min_res = std::min(fabs(res), min_res);

                if (fabs(res) > EPSILON)
                    out_plan.push_back(face_vertexes[i]);
                else
                    inside_plan.push_back(face_vertexes[i]);
            }
            if (out_plan.size() > 0)
                std::cout << "shape " << num_shp << ", face " << f << ", inside " << inside_plan.size() << ", outside " << out_plan.size() << ", max " << max_res << ", min " << min_res << std::endl;
        }
}

void tessToshp::massProperties() {
    for (size_t i = 0; i < Shapes.size(); i++) {
        Shapes[i]->defineVertexConnectivity();
        Shapes[i]->MCnstep = 50000;
        Shapes[i]->massProperties();
        Shapes[i]->preCompDone = 'y';
    }

    for (size_t i = 0; i < Shapes.size(); i++) {
        Shapes[i]->fitObb();
        int iv_out = -1;
        do {
            if (iv_out > -1) {
                Shapes[i]->obb.enlarge(radius);
                iv_out = -1;
            }
        } while (!all_vertexes_in_shape(*Shapes[i], iv_out));
    }
}

bool tessToshp::all_vertexes_in_shape(Shape &shape, int &iv_out) {
    OBB obb = shape.obb;
    //obb.enlarge(-shape.radius);
    for (size_t iv = 0; iv < shape.vertex.size(); iv++) {
        vec3r v = shape.vertex[iv];
        if (!obb.intersect(v)) {
            iv_out = (int)iv;
            return false;
        }
    }
    return true;
}

void tessToshp::writeShapes() {
    double maxRadius = -1.e+30;
    double minRadius = +1.e+30;

    dVerlet = +1.e+30;
    DVerlet = -1.e+30;

    std::string fout_ = "./shape.shp";
    std::ofstream out(fout_.c_str());

    shapeFile = fout_.c_str();

    for (size_t s = 0; s < Shapes.size(); ++s) {
        Shape &shp = *Shapes[s];
        if (shp.obb.extent.x == shp.obb.extent.y && shp.obb.extent.y == shp.obb.extent.z && shp.obb.extent.z == radius) {
            OBB &obb = shp.obb;
            AABB aabb;
            shp.getAABB(aabb);

            obb.e[0] = vec3r::unit_x();
            obb.e[1] = vec3r::unit_y();
            obb.e[2] = vec3r::unit_z();
            obb.extent = (aabb.max - aabb.min) / 2.;
            obb.center = (aabb.max + aabb.min) / 2.;
        }
    }

    for (size_t s = 0; s < Shapes.size(); ++s) {
        dVerlet = std::min(dVerlet, Shapes[s]->radius);
        DVerlet = std::max(DVerlet, Shapes[s]->radius);
    }

    for (size_t s = 0; s < Shapes.size(); ++s) {
        minRadius = std::min(minRadius, Shapes[s]->obb.extent.length());
        maxRadius = std::max(maxRadius, Shapes[s]->obb.extent.length());
    }

    std::cout << "maxRadius " << maxRadius << " minRadius " << minRadius << std::endl;

    dVerlet = std::min(dVerlet, minRadius);
    DVerlet = std::max(DVerlet, maxRadius);

    std::cout << "DVerlet " << DVerlet << " dVerlet " << dVerlet << std::endl;

    Particles.resize(Shapes.size());

    for (size_t i = 0; i < Particles.size(); i++) {
        Particle &P = Particles[i];
        Particles[i].shape = (Shapes[i]);
        quat Q = P.shape->orientation;
        P.Q = Q;
        P.pos = P.shape->position;
        P.homothety = 1.;
    }

    if (out) {
        for (size_t s = 0; s < Shapes.size(); ++s) {
            CalculateAreas(*Shapes[s]);
            Shapes[s]->write(out);
        }
        out.close();
    } else {
        std::cerr << "Error message : "
                  << "Impossible de lire le fichier : " << fout_.c_str() << std::endl;
    }
}

void tessToshp::CalculateAreas(Shape &shp) {
    shp.faceArea.resize(shp.face.size());
    shp.weight.resize(shp.face.size());

    for (size_t f = 0; f < shp.face.size(); f++) {
        //centre O
        vec3r O = vec3r();
        size_t N = shp.face[f].size();
        shp.weight[f].resize(N);

        for (size_t iv = 0; iv < shp.face[f].size(); iv++) {
            vec3r v = shp.vertex[shp.face[f][iv]];
            O += v;
        }

        O /= N;

        double area = 0;
        vec3r Area;

        for (size_t iv = 0; iv < shp.face[f].size(); iv++) {
            vec3r A = shp.vertex[shp.face[f][iv]];
            vec3r B = shp.vertex[shp.face[f][(iv + 1) % N]];
            vec3r OA = A - O;
            vec3r OB = B - O;
            Area += .5 * (OA ^ OB);
        }

        area = norm(Area);
        shp.faceArea[f] = area;

        shp.weight[f].resize(shp.face[f].size());

        for (size_t iv = 0; iv < shp.face[f].size(); iv++) {
            /* M-------A-------P */

            int ivm = iv - 1;
            if (ivm < 0)
                ivm = N - 1;
            int ivp = iv + 1;
            if (ivp > (int)(N - 1))
                ivp = 0;

            vec3r Am = shp.vertex[shp.face[f][ivm]];
            vec3r A = shp.vertex[shp.face[f][iv]];
            vec3r Ap = shp.vertex[shp.face[f][ivp]];

            vec3r M = (A + Am) / 2.;
            vec3r P = (A + Ap) / 2.;

            vec3r OM = M - O;
            vec3r OA = A - O;
            vec3r OP = P - O;

            vec3r w = (OM ^ OA) + (OA ^ OP);

            if (area != 0.)
                shp.weight[f][iv] = 0.5 * norm(w) / area;
            else
                shp.weight[f][iv] = 0.;
        }
    }
}

/*
const Vec2 findCentroid(Vec2* pts, size_t nPts){
    Vec2 off = pts[0];
    float twicearea = 0;
    float x = 0;
    float y = 0;
    Vec2 p1, p2;
    float f;
    for (int i = 0, j = nPts - 1; i < nPts; j = i++) {
        p1 = pts[i];
        p2 = pts[j];
        f = (p1.x - off.x) * (p2.y - off.y) - (p2.x - off.x) * (p1.y - off.y);
        twicearea += f;
        x += (p1.x + p2.x - 2 * off.x) * f;
        y += (p1.y + p2.y - 2 * off.y) * f;
    }

    f = twicearea * 3;

    return Vec2(x / f + off.x, y / f + off.y);
}
*/

/*
import numpy as np
#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
x = np.linalg.det([[1,a[1],a[2]],
[1,b[1],b[2]],
[1,c[1],c[2]]])
y = np.linalg.det([[a[0],1,a[2]],
[b[0],1,b[2]],
[c[0],1,c[2]]])
z = np.linalg.det([[a[0],a[1],1],
[b[0],b[1],1],
[c[0],c[1],1]])
magnitude = (x**2 + y**2 + z**2)**.5
return (x/magnitude, y/magnitude, z/magnitude)
#area of polygon poly
def poly_area(poly):
if len(poly) < 3: # not a plane - no area
return 0
total = [0, 0, 0]
N = len(poly)
for i in range(N):
vi1 = poly[i]
vi2 = poly[(i+1) % N]
prod = np.cross(vi1, vi2)
total[0] += prod[0]
total[1] += prod[1]
total[2] += prod[2]
result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
return abs(result/2)
*/

void tessToshp::loadConfNeper(const char *name) {
    std::ifstream conf(name);
    if (!conf.is_open()) {
        std::cerr << msg::warn() << "@Rockable, Cannot read " << name << msg::normal() << std::endl;
        return;
    }

    // Check header
    std::string prog;
    conf >> prog;
    if (prog != "Rockable") {
        std::cerr << msg::warn() << "@Rockable, This is not a file for the code Rockable!" << msg::normal() << std::endl;
    }
    std::string date;
    conf >> date;
    if (date != CONF_VERSION_DATE) {
        std::cerr << msg::warn() << "@Rockable, The version-date should be " << CONF_VERSION_DATE << "!" << msg::normal()
                  << std::endl;
    }

    kwParser parser;

    parser.kwMap["t"] = __GET__(conf, t);
    parser.kwMap["tmax"] = __GET__(conf, tmax);
    parser.kwMap["dt"] = __DO__(conf) {
        conf >> dt;
        dt_2 = 0.5 * dt;
        dt2_2 = 0.5 * dt * dt;
    };
    parser.kwMap["interVerlet"] = __GET__(conf, interVerlet);
    parser.kwMap["interConf"] = __GET__(conf, interConf);

    parser.kwMap["forceLaw"] = __DO__(conf) {
        std::string lawName;
        conf >> lawName;
        if (lawName == "Avalanches") {
            forceLawPtr = bind(&Rockable::forceLawAvalanches, this, std::placeholders::_1);
            optionNames["forceLaw"] = "Avalanches";
        } else if (lawName == "StickedLinks") {
            forceLawPtr = bind(&Rockable::forceLawStickedLinks, this, std::placeholders::_1);
            optionNames["forceLaw"] = "StickedLinks";
        } else
            std::cout << msg::warn() << "forceLaw " << lawName << " is unknown" << msg::normal() << std::endl;
    };

  parser.kwMap["knContact"] = __DO__(conf) { readLawData(conf, idKnContact); };
  parser.kwMap["en2Contact"] = __DO__(conf) { readLawData(conf, idEn2Contact); };
  parser.kwMap["ktContact"] = __DO__(conf) { readLawData(conf, idKtContact); };
  parser.kwMap["muContact"] = __DO__(conf) { readLawData(conf, idMuContact); };
  parser.kwMap["krContact"] = __DO__(conf) { readLawData(conf, idKrContact); };
  parser.kwMap["murContact"] = __DO__(conf) { readLawData(conf, idMurContact); };

    parser.kwMap["iconf"] = __GET__(conf, iconf);
    parser.kwMap["nDriven"] = __GET__(conf, nDriven);
    parser.kwMap["Sphere_Radius"] = __GET__(conf, R);
    parser.kwMap["Cuts_Number"] = __GET__(conf, nbP);
    parser.kwMap["Cells_Number"] = __GET__(conf, number_of_cells);
    parser.kwMap["MinskowskiRadius"] = __GET__(conf, radius);

    // This single line actually parses the file
    parser.parse(conf);
    if (conf.is_open())
        conf.close();
}

//LHASSAN
/**
 @brief Save a input.txt '
 */
void tessToshp::saveConfNeper() {
    char fname[256];
    sprintf(fname, "input.txt");
    std::ofstream conf(fname);

    conf << "Rockable " << CONF_VERSION_DATE << '\n'; // format: progName version-date(dd-mm-yyyy)
    conf << "t " << t << '\n';
    conf << "tmax " << tmax << '\n';
    conf << "dt " << dt << '\n';
    conf << "interVerlet " << interVerlet << '\n';
    conf << "interConf " << interConf << '\n';
    conf << "DVerlet " << DVerlet << '\n';
    conf << "dVerlet " << dVerlet << '\n';
    for (size_t grp = 0; grp < properties.ngroup; grp++) {
        double density = properties.get(idDensity, grp);
        if (density > 0.0) // TODO implement properties.isDefined(idDensity, grp)
            conf << "density " << grp << " " << density << '\n';
    }
    conf << "gravity " << gravity << '\n';
    if (bodyForce != nullptr) {
        conf << "BodyForce ";
        bodyForce->write(conf);
    }

    for (auto it = optionNames.begin(); it != optionNames.end(); ++it) {
        if (it->first == "forceLaw")
            conf << it->first << " " << it->second << '\n';
    }

    writeLawData(conf, "knContact");
    writeLawData(conf, "en2Contact");
    writeLawData(conf, "ktContact");
    writeLawData(conf, "muContact");
    writeLawData(conf, "krContact");
    writeLawData(conf, "murContact");

    conf << "iconf " << iconf << '\n';
    conf << "nDriven " << nDriven << '\n';
    conf << "shapeFile " << shapeFile << '\n';
    if (CommBox().sep != ' ')
        conf << "separator " << CommBox().keywordFromSep() << '\n';
    conf << std::scientific << std::setprecision(CommBox().precision);
    conf << "Particles " << Particles.size() << '\n';
    conf << "#name" << CommBox().sep << "group" << CommBox().sep << "cluster" << CommBox().sep << "homothety"
         << CommBox().sep << "pos.x" << CommBox().sep << "pos.y" << CommBox().sep << "pos.z" << CommBox().sep << "vel.x"
         << CommBox().sep << "vel.y" << CommBox().sep << "vel.z" << CommBox().sep << "acc.x" << CommBox().sep << "acc.y"
         << CommBox().sep << "acc.z" << CommBox().sep << "Q.w" << CommBox().sep << "Q.x" << CommBox().sep << "Q.y"
         << CommBox().sep << "Q.z" << CommBox().sep << "vrot.x" << CommBox().sep << "vrot.y" << CommBox().sep << "vrot.z"
         << CommBox().sep << "arot.x" << CommBox().sep << "arot.y" << CommBox().sep << "arot.z" << '\n';

    for (size_t i = 0; i < Particles.size(); i++) {
        conf << Particles[i].shape->name << CommBox().sep << Particles[i].group << CommBox().sep << Particles[i].cluster
             << CommBox().sep << Particles[i].homothety << CommBox().sep << Particles[i].pos << CommBox().sep
             << Particles[i].vel << CommBox().sep << Particles[i].acc << CommBox().sep << Particles[i].Q << CommBox().sep
             << Particles[i].vrot << CommBox().sep << Particles[i].arot << '\n';
    }
    conf << std::flush;

    if (conf.is_open())
        conf.close();
}

void tessToshp::exec() {
    min_norm_edge = +1.e+20;
    max_norm_edge = -1.e+20;
    //createShpere();
    createTess();
    readTess();
    extractShapes();
    for (size_t s = 0; s < Shapes.size(); s++)
        validate_shape_faces(s);

    std::cout << "_________________ Before _________________" << std::endl;
    std::cout << " MIN NROM EDGE " << min_norm_edge << std::endl;
    std::cout << " MAX NROM EDGE " << max_norm_edge << std::endl;
    std::cout << "______________________________________" << std::endl;

    create_shapes();

    std::cout << "_________________ After _________________" << std::endl;
    std::cout << " MIN NROM EDGE " << min_norm_edge << std::endl;
    std::cout << " MAX NROM EDGE " << max_norm_edge << std::endl;
    std::cout << "______________________________________" << std::endl;

    massProperties();
    writeShapes();
}
