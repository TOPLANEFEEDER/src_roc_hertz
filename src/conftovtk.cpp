//  Copyright or © or Copr. Rockable
//
//  vincent.richefeu@3sr-grenoble.fr
//
//  This software is a computer program whose purpose is
//    (i)  to hold sphero-polyhedral shapes,
//    (ii) to manage breakable interfaces.
//  It is developed for an ACADEMIC USAGE
//
//  This software is governed by the CeCILL-B license under French law and
//  abiding by the rules of distribution of free software.  You can  use,
//  modify and/ or redistribute the software under the terms of the CeCILL-B
//  license as circulated by CEA, CNRS and INRIA at the following URL
//  "http://www.cecill.info".
//
//  As a counterpart to the access to the source code and  rights to copy,
//  modify and redistribute granted by the license, users are provided only
//  with a limited warranty  and the software's author,  the holder of the
//  economic rights,  and the successive licensors  have only  limited
//  liability.
//
//  In this respect, the user's attention is drawn to the risks associated
//  with loading,  using,  modifying and/or developing or reproducing the
//  software by the user in light of its specific status of free software,
//  that may mean  that it is complicated to manipulate,  and  that  also
//  therefore means  that it is reserved for developers  and  experienced
//  professionals having in-depth computer knowledge. Users are therefore
//  encouraged to load and test the software's suitability as regards their
//  requirements in conditions enabling the security of their systems and/or
//  data to be ensured and,  more generally, to use and operate it in the
//  same conditions as regards security.
//
//  The fact that you are presently reading this means that you have had
//  knowledge of the CeCILL-B license and that you accept its terms.

#include "conftovtk.hpp"

bool tryToReadConf(int num)
{
    char file_name[256];
    sprintf(file_name, "conf%d", num);
    if (fileTool::fileExists(file_name))
    {
        std::cout << "Read " << file_name << std::endl;
        box.clearMemory();
        box.loadConf(file_name);
        complexityNumber = 0;
        for (size_t i = 0; i < box.Particles.size(); ++i)
        {
            complexityNumber += box.Particles[i].shape->vertex.size();
        }
        confNum = box.iconf;
        box.computeAABB();
    }
    else
    {
        std::cout << file_name << " does not exist" << std::endl;
        return false;
    }
    return true;
}

bool iscolinearity(const vec3r &v0, const vec3r &v1, const vec3r &v2)
{
    vec3r v01 = v1 - v0;
    v01.normalized();
    vec3r v02 = v2 - v0;
    v02.normalized();
    double cos0 = dot(v01, v02);
    double eps = 1.e-10;
    return ((1. - fabs(cos0)) < eps);
}

std::string TruncateLongString(std::string const str, unsigned int maxLength)
{
    std::string str_return = str;

    if (str.length() > maxLength)
        return str_return.substr(0, maxLength);
    return str_return;
}

void SeparateParticlesByType()
{
    Sphers.clear();
    Joncs.clear();
    Polyrs.clear();

    SphersId.clear();
    JoncsId.clear();
    PolyrsId.clear();

    for (size_t id = 0; id < box.Particles.size(); id++)
    {
        Particle *P = &box.Particles[id];

        if (P->shape->vertex.size() == 1)
        {
            Sphers.push_back(P);
            SphersId.push_back(id);
        }
        if (P->shape->vertex.size() == 2)
        {
            Joncs.push_back(P);
            JoncsId.push_back(id);
        }
        if (P->shape->vertex.size() > 2)
        {
            Polyrs.push_back(P);
            PolyrsId.push_back(id);
        }
    }
}

void writeVTKSPHER(int num)
{
    unsigned int i;
    char vtk_file[200];
    FILE *sortie_vtk;

    unsigned int nbgrains = (unsigned int)Sphers.size();

    sprintf(vtk_file, "Shepres%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");

    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "Sortie Particles\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(sortie_vtk, "POINTS %i float\n", nbgrains);

    for (i = 0; i < nbgrains; i++)
    {
        vec3r &X = Sphers[i]->pos;
        fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }

    fprintf(sortie_vtk, "POINT_DATA %i\n", nbgrains);

    fprintf(sortie_vtk, "VECTORS Radius float\n");

    for (i = 0; i < nbgrains; i++)
    {
        fprintf(sortie_vtk, "0. 0. %lf\n", Sphers[i]->MinskowskiRadius());
    }

    fprintf(sortie_vtk, "VECTORS Vel float\n");

    for (i = 0; i < nbgrains; i++)
    {
        vec3r &V = Sphers[i]->vel;
        fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
    }

    fprintf(sortie_vtk, "SCALARS Material float\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

    for (i = 0; i < nbgrains; i++)
    {
        fprintf(sortie_vtk, "%i\n", Sphers[i]->group);
    }

    fprintf(sortie_vtk, "SCALARS Id int\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

    for (i = 0; i < nbgrains; i++)
    {
        fprintf(sortie_vtk, "%li\n", SphersId[i]);
    }

    fclose(sortie_vtk);
}

/*
 // Draw the shape of the sphero-polyhedron in its own framework
void drawShape(Shape* s, double homothety) {
  double R = homothety * s->radius;
  int nbLevelSphere = 2;
  if (complexityNumber > 10000) {
    nbLevelSphere = 1;
  }

  if (complexMode == 0) {
    // vertixes (spheres)
    for (size_t v = 0; v < s->vertex.size(); ++v) {
      glPushMatrix();
      glTranslatef(homothety * s->vertex[v].x, homothety * s->vertex[v].y, homothety * s->vertex[v].z);
      facetSphere::draw(nbLevelSphere, R);
      glPopMatrix();
    }

    // edges (tubes)
    for (size_t e = 0; e < s->edge.size(); ++e) {
      vec3r orig = homothety * s->vertex[s->edge[e].first];
      vec3r arrow = homothety * s->vertex[s->edge[e].second];
      arrow -= orig;
      glutShape::drawTube(orig, arrow, 2.0 * R);
    }
  }

  // faces (3D polygones)
  for (size_t f = 0; f < s->face.size(); ++f) {
    if (s->face[f].size() < 3) continue;  // At least 3 pts!
    vec3r N =
        cross(s->vertex[s->face[f][1]] - s->vertex[s->face[f][0]], s->vertex[s->face[f][2]] - s->vertex[s->face[f][0]]);
    N.normalize();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(N.x, N.y, N.z);
    for (size_t v = 0; v < s->face[f].size(); ++v) {
      glVertex3f(homothety * s->vertex[s->face[f][v]].x + N.x * R, homothety * s->vertex[s->face[f][v]].y + N.y * R,
                 homothety * s->vertex[s->face[f][v]].z + N.z * R);
    }
    glEnd();

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(-N.x, -N.y, -N.z);
    for (size_t v = 0; v < s->face[f].size(); ++v) {
      glVertex3f(homothety * s->vertex[s->face[f][v]].x - N.x * R, homothety * s->vertex[s->face[f][v]].y - N.y * R,
                 homothety * s->vertex[s->face[f][v]].z - N.z * R);
    }
    glEnd();
  }
}

void drawParticles() {
  if (mouse_mode != NOTHING && complexityNumber > 5000) {
    drawOBBs();
    return;
  }

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  for (size_t i = box.nDriven; i < box.Particles.size(); ++i) {
    if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
      glColor4ub(234, 255, 0, (int)floor(alpha_particles * 255));
    } else
      glColor4ub(178, 34, 34, (int)floor(alpha_particles * 255));

    vec3r pos = box.Particles[i].pos;

    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
    glMultMatrixf(Rot_Matrix);
    drawShape(box.Particles[i].shape, box.Particles[i].homothety);
    glPopMatrix();
  }

  for (size_t i = 0; i < box.nDriven; ++i) {
    if (selectedParticle >= 0 && i == (size_t)selectedParticle) {
      glColor4ub(234, 255, 0, (int)floor(alpha_fixparticles * 255));
    } else
      glColor4ub(128, 128, 128, (int)floor(alpha_fixparticles * 255));

    vec3r pos = box.Particles[i].pos;

    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    quat2GLMatrix<GLfloat>(box.Particles[i].Q, Rot_Matrix);
    glMultMatrixf(Rot_Matrix);
    drawShape(box.Particles[i].shape, box.Particles[i].homothety);
    glPopMatrix();
  }
}
*/

void writeVTKPOLYR(int num)
{
    char vtk_file[200];
    FILE *sortie_vtk;
    int vertex_compt = 0;
    int nb_points = 0;
    int nb_faces = 0;
    int nb_face_points = 0;
    unsigned int nbgrains = (unsigned int)Polyrs.size();
    for (size_t i = 0; i < nbgrains; i++)
    {
        nb_points += Polyrs[i]->shape->vertex.size();
        nb_faces += Polyrs[i]->shape->face.size();

        for (size_t j = 0; j < Polyrs[i]->shape->face.size(); j++)
        {
            nb_face_points += Polyrs[i]->shape->face[j].size();
        }
    }

    sprintf(vtk_file, "Polyr%.4i.vtk", num);

    sortie_vtk = fopen(vtk_file, "w");

    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "RIGID      1\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET POLYDATA\n");
    fprintf(sortie_vtk, "POINTS %i float\n", nb_points);

    // # Ecriture des Vertices
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;

        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            vec3r X = Polyrs[i]->Glob(shp->vertex[j]);
            vec3r n = (X - Polyrs[i]->pos);
            n.normalize();
            X += n * shp->radius * Polyrs[i]->homothety;
            fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
    }

    // Ecriture des Faces (connectivity des vertices)
    fprintf(sortie_vtk, "POLYGONS %i %i\n", nb_faces, nb_face_points + nb_faces);

    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%li ", shp->face[j].size());
            for (size_t k = 0; k < shp->face[j].size(); k++)
                fprintf(sortie_vtk, "%li ", vertex_compt + shp->face[j][k]);
            fprintf(sortie_vtk, "\n");
        }
        vertex_compt += shp->vertex.size();
    }

    fprintf(sortie_vtk, "CELL_DATA %i \n", nb_faces);
    fprintf(sortie_vtk, "SCALARS Material int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%d\n", Polyrs[i]->group);
        }
    }

    fprintf(sortie_vtk, "SCALARS Id int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%li\n", PolyrsId[i]);
        }
    }


    fprintf(sortie_vtk, "VECTORS Vel float\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        vec3r &V = Polyrs[i]->vel;
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
        }
    }

    fclose(sortie_vtk);
}

//@dv:Elongated particles
void writeVTKJONCx(int num)
{
    char vtk_file[200];
    char vtk_file2[200];
    FILE *sortie_vtk;
    FILE *sortie_vtk2;
    int vertex_compt = 0;
    int nb_points = 0;
    int nb_edges = 0;
    unsigned int nbgrains = (unsigned int)Joncs.size();
    for (size_t i = 0; i < nbgrains; i++)
    {
        nb_points += Joncs[i]->shape->vertex.size();
        nb_edges += Joncs[i]->shape->edge.size();
    }

    sprintf(vtk_file, "Jonc%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");
    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "RIGID      1\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET POLYDATA\n");
    fprintf(sortie_vtk, "POINTS %i float\n", nb_points);

    sprintf(vtk_file2, "JoncSpheres%.4i.vtk", num);
    sortie_vtk2 = fopen(vtk_file2, "w");
    fprintf(sortie_vtk2, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk2, "Sortie Particles\n");
    fprintf(sortie_vtk2, "ASCII\n");
    fprintf(sortie_vtk2, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(sortie_vtk2, "POINTS %i float\n", nb_points);

    // # Ecriture des Vertices
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;

        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            vec3r X = Joncs[i]->Glob(shp->vertex[j]);
            fprintf(sortie_vtk2, "%lf %lf %lf\n", X.x, X.y, X.z);
            // vec3r n = (X - Joncs[i]->pos);
            // n.normalize();
            // X += n*shp->radius*Joncs[i]->homothety;
            fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
    }

    fprintf(sortie_vtk2, "POINT_DATA %i\n", nb_points);
    fprintf(sortie_vtk2, "VECTORS Radius float\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;

        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            fprintf(sortie_vtk2, "0. 0. %lf\n", shp->radius);
        }
    }

    fprintf(sortie_vtk2, "SCALARS Id int 1\n");
    fprintf(sortie_vtk2, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;
        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            fprintf(sortie_vtk2, "%li\n", JoncsId[i]);
        }
    }

    fclose(sortie_vtk2);

    // Ecriture des Edges (connectivity des vertices)

    fprintf(sortie_vtk, "LINES %i %i\n", nb_edges, 3 * nb_edges);

    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;
        size_t nbe = shp->edge.size();
        for (size_t j = 0; j < nbe; j++)
        {
            fprintf(sortie_vtk, "2 %li %li\n", vertex_compt + shp->edge[j].first, vertex_compt + shp->edge[j].second);
        }
        vertex_compt += shp->vertex.size();
    }

    fprintf(sortie_vtk, "CELL_DATA %i \n", nb_edges);
    fprintf(sortie_vtk, "SCALARS Material int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;
        size_t nbe = shp->edge.size();
        for (size_t j = 0; j < nbe; j++)
        {
            fprintf(sortie_vtk, "%d\n", Joncs[i]->group);
        }
    }

    fprintf(sortie_vtk, "SCALARS Id int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Joncs[i]->shape;
        size_t nbe = shp->edge.size();
        for (size_t j = 0; j < nbe; j++)
        {
            fprintf(sortie_vtk, "%li\n", JoncsId[i]);
        }
    }

    fprintf(sortie_vtk, "VECTORS Vel float\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        vec3r &V = Joncs[i]->vel;
        Shape *shp = Joncs[i]->shape;
        size_t nbe = shp->edge.size();
        for (size_t j = 0; j < nbe; j++)
        {
            fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
        }
    }

    fclose(sortie_vtk);
}

void writeVTKPOLYRVertex(int num)
{
    int vertex_compt = 0;
    int nb_points = 0;
    int nb_faces = 0;
    int nb_face_points = 0;
    unsigned int nbgrains = (unsigned int)Polyrs.size();
    for (size_t i = 0; i < nbgrains; i++)
    {
        nb_points += Polyrs[i]->shape->vertex.size();
        nb_faces += Polyrs[i]->shape->face.size();

        for (size_t j = 0; j < Polyrs[i]->shape->face.size(); j++)
        {
            nb_face_points += Polyrs[i]->shape->face[j].size();
        }
    }

    char vtk_file[200];
    char vtk_file2[200];

    FILE *sortie_vtk;
    FILE *sortie_vtk2;

    sprintf(vtk_file, "Polyr%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");
    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "RIGID      1\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET POLYDATA\n");
    fprintf(sortie_vtk, "POINTS %i float\n", nb_points);

    sprintf(vtk_file2, "PolyrSpheres%.4i.vtk", num);
    sortie_vtk2 = fopen(vtk_file2, "w");
    fprintf(sortie_vtk2, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk2, "Sortie Particles\n");
    fprintf(sortie_vtk2, "ASCII\n");
    fprintf(sortie_vtk2, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(sortie_vtk2, "POINTS %i float\n", nb_points);

    // # Ecriture des Vertices
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;

        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            vec3r X = Polyrs[i]->Glob(shp->vertex[j]);
            fprintf(sortie_vtk2, "%lf %lf %lf\n", X.x, X.y, X.z);
            // vec3r n = (X - Polyrs[i]->pos);
            // n.normalize();
            // X += 2.*n*shp->radius*Polyrs[i]->homothety;
            fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
        }
    }

    fprintf(sortie_vtk2, "POINT_DATA %i\n", nb_points);
    fprintf(sortie_vtk2, "VECTORS Radius float\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;

        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            fprintf(sortie_vtk2, "0. 0. %lf\n", shp->radius);
        }
    }

    fprintf(sortie_vtk2, "SCALARS Id int 1\n");
    fprintf(sortie_vtk2, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbv = shp->vertex.size();
        for (size_t j = 0; j < nbv; j++)
        {
            fprintf(sortie_vtk2, "%li\n", PolyrsId[i]);
        }
    }

    fclose(sortie_vtk2);

    // Ecriture des Faces (connectivity des vertices)
    fprintf(sortie_vtk, "POLYGONS %i %i\n", nb_faces, nb_face_points + nb_faces);

    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%li ", shp->face[j].size());
            for (size_t k = 0; k < shp->face[j].size(); k++)
                fprintf(sortie_vtk, "%li ", vertex_compt + shp->face[j][k]);
            fprintf(sortie_vtk, "\n");
        }
        vertex_compt += shp->vertex.size();
    }

    fprintf(sortie_vtk, "CELL_DATA %i \n", nb_faces);
    fprintf(sortie_vtk, "SCALARS Material int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%d\n", Polyrs[i]->group);
        }
    }

    fprintf(sortie_vtk, "SCALARS Id int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%li\n", PolyrsId[i]);
        }
    }

    fprintf(sortie_vtk, "VECTORS Vel float\n");
    for (size_t i = 0; i < nbgrains; i++)
    {
        vec3r &V = Polyrs[i]->vel;
        Shape *shp = Polyrs[i]->shape;
        size_t nbf = shp->face.size();
        for (size_t j = 0; j < nbf; j++)
        {
            fprintf(sortie_vtk, "%lf %lf %lf\n", V.x, V.y, V.z);
        }
    }

    fclose(sortie_vtk);
}

void writeVTKContactsShperes(int num)
{
    char vtk_file[200];
    FILE *sortie_vtk;

    sprintf(vtk_file, "Forces%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");

    long nbContacts = 0;
    double fMax = 0;

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    nbContacts++;
                    fMax = std::max(fMax, norm(it->fn * it->n + it->ft));
                }
            }
        }
    }

    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "Sortie Particles\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(sortie_vtk, "POINTS %li float\n", nbContacts);

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    vec3r X = it->pos;
                    fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
                }
            }
        }
    }

    fprintf(sortie_vtk, "POINT_DATA %li\n", nbContacts);

    fprintf(sortie_vtk, "VECTORS F_Fmax float\n");

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                int i = it->i;
                int j = it->j;

                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    fMax = box.Particles[i].MinskowskiRadius() + box.Particles[j].MinskowskiRadius();
                    if (fMax != 0.)
                        fprintf(sortie_vtk, "0. 0. %lf\n", fMax / 2.);
                    else
                        fprintf(sortie_vtk, "0. 0. 0.\n");
                }
            }
        }
    }

    fprintf(sortie_vtk, "VECTORS Fn float\n");

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    vec3r X = it->fn * it->n;
                    fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
                }
            }
        }
    }

    fprintf(sortie_vtk, "VECTORS Ft float\n");

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    vec3r X = it->ft;
                    fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
                }
            }
        }
    }

    fprintf(sortie_vtk, "SCALARS Type float\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        std::set<Interaction>::iterator it = box.Interactions[k].begin();
        for (; it != box.Interactions[k].end(); ++it)
        {
            {
                if (it->dn > 0.0 && it->stick == nullptr)
                    continue;
                else
                {
                    fprintf(sortie_vtk, "%d\n", it->type);
                }
            }
        }
    }

    fclose(sortie_vtk);
}

//@dv
void writeVTKContactsPolys(int num)
{
    char vtk_file[200];
    FILE *sortie_vtk;

    sprintf(vtk_file, "Forces%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");

    long nbContacts = 0;
    double fMax = 0.0;

    //@dv
    box.activeInteractions.clear();

    for (size_t k = 0; k < box.Interactions.size(); ++k)
    {
        for (auto it = box.Interactions[k].begin(); it != box.Interactions[k].end(); ++it)
        {
            Interaction *I = const_cast<Interaction *>(std::addressof(*it));
            if (it->dn < 0.0)
            {
                box.activeInteractions.push_back(I);
            }
        }
    }

    for (size_t k = 0; k < box.activeInteractions.size(); ++k)
    {
        Interaction *I = box.activeInteractions[k];
        if (I->i > I->j)
        {
            Interaction other = *I;

            I->i = other.j;
            I->j = other.i;
            I->type = other.type;
            I->isub = other.jsub;
            I->jsub = other.isub;
            I->prev_n = -other.prev_n;
            I->n = -other.n;
            I->dn = other.dn;
            I->prev_dn = other.prev_dn;
            I->pos = other.pos;
            I->posI = other.posJ; //@dv:change I->J
            I->posJ = other.posI;
            I->lji = -other.lji;
            I->vel = -other.vel;
            I->fn = other.fn;
            I->ft = -other.ft;
            I->mom = -other.mom;
            I->damp = other.damp;
            I->stick = other.stick;
        }
    }

    std::sort(box.activeInteractions.begin(), box.activeInteractions.end(), std::less<Interaction *>());

    size_t ncont = 0;

    size_t i = box.activeInteractions[0]->i;
    size_t j = box.activeInteractions[0]->j;

    size_t ibefore = i;
    size_t jbefore = j;

    contacts.clear();
    contacts.push_back(new Contact);
    contacts[ncont]->i = i;
    contacts[ncont]->j = j;

    for (size_t k = 0; k < box.activeInteractions.size(); ++k)
    {
        Interaction *I = box.activeInteractions[k];

        i = I->i;
        j = I->j;

        if (i != ibefore || j != jbefore) // i or j changes, so we add new Contact to the contacts vector
        {
            contacts.push_back(new Contact);
            ncont++; // we increment ncont after push_back
            contacts[ncont]->i = i;
            contacts[ncont]->j = j;
            contacts[ncont]->subinter.push_back(I); // we add the interaction to the subinter
            ibefore = i;
            jbefore = j;
        }
        else // i and j doesn't change
        {
            contacts[ncont]->subinter.push_back(I); // here i and j dont changes so we just add the interaction to the subinter
        }
    }

    // we calculate the fn and ft vector (vectorial sum)
    for (size_t c = 0; c < contacts.size(); ++c)
    {
        Contact *cij = contacts[c];

        cij->fn.reset();
        cij->ft.reset();

        for (size_t k = 0; k < cij->subinter.size(); ++k)
        {
            Interaction *I = cij->subinter[k];
            cij->fn += I->fn * I->n;
            cij->ft += I->ft;
        }

        // fMax = std::max(fMax, norm(cij->fn + cij->ft));
    }

    nbContacts = contacts.size();

    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "Sortie Particles\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(sortie_vtk, "POINTS %li float\n", nbContacts);

    for (size_t c = 0; c < contacts.size(); c++)
    {
        Contact *cij = contacts[c];
        cij->posI.reset();
        cij->posJ.reset();
        for (size_t k = 0; k < cij->subinter.size(); k++)
        {
            Interaction *I = cij->subinter[k];
            cij->posI += I->posI;
            cij->posJ += I->posJ;
        }
        cij->posI /= cij->subinter.size();
        cij->posJ /= cij->subinter.size();

        vec3r X = cij->posI;
        fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }

    fprintf(sortie_vtk, "POINT_DATA %li\n", nbContacts);
    fprintf(sortie_vtk, "VECTORS F_Fmax float\n");

    for (size_t c = 0; c < contacts.size(); c++)
    {
        Contact *cij = contacts[c];
        cij->l = cij->subinter[0]->lji;
        fMax = norm(cij->l);
        if (fMax != 0.)
            fprintf(sortie_vtk, "0. 0. %lf\n", fMax / 2.);
        else
            fprintf(sortie_vtk, "0. 0. 0.\n");
    }

    fprintf(sortie_vtk, "VECTORS Fn float\n");

    for (size_t c = 0; c < contacts.size(); c++)
    {
        Contact *cij = contacts[c];
        vec3r X = cij->fn;
        fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }

    fprintf(sortie_vtk, "VECTORS Ft float\n");

    for (size_t c = 0; c < contacts.size(); c++)
    {
        Contact *cij = contacts[c];
        vec3r X = cij->ft;
        fprintf(sortie_vtk, "%lf %lf %lf\n", X.x, X.y, X.z);
    }

    fclose(sortie_vtk);
}


void writeVTKOBB(int num)
{
    char vtk_file[200];
    FILE *sortie_vtk;

    sprintf(vtk_file, "OBB%.4i.vtk", num);
    sortie_vtk = fopen(vtk_file, "w");

    size_t nbSommets = 8 * box.Particles.size();

    fprintf(sortie_vtk, "# vtk DataFile Version 3.0\n");
    fprintf(sortie_vtk, "Sortie Foces\n");
    fprintf(sortie_vtk, "ASCII\n");
    fprintf(sortie_vtk, "DATASET POLYDATA\n");

    fprintf(sortie_vtk, "POINTS %li float\n", nbSommets);

    OBB obbi;

    for (size_t i = 0; i < box.Particles.size(); i++)
    {

        obbi = box.Particles[i].shape->obb;
        obbi.rotate(box.Particles[i].Q);
        obbi.extent *= box.Particles[i].homothety;
        obbi.center *= box.Particles[i].homothety;
        obbi.center += box.Particles[i].pos;
        // if (enlarged_obb) obbi.enlarge(0.5 * box.DVerlet);

        // je vais écrire les 8 sommets
        vec3r corner;

        // les 4 sommets du bas
        corner = obbi.center - obbi.extent[0] * obbi.e[0] - obbi.extent[1] * obbi.e[1] - obbi.extent[2] * obbi.e[2];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner += 2.0 * obbi.extent[0] * obbi.e[0];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner += 2.0 * obbi.extent[1] * obbi.e[1];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner -= 2.0 * obbi.extent[0] * obbi.e[0];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);

        // les 4 sommets du haut
        corner = obbi.center - obbi.extent[0] * obbi.e[0] - obbi.extent[1] * obbi.e[1] - obbi.extent[2] * obbi.e[2];
        corner += 2.0 * obbi.extent[2] * obbi.e[2];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner += 2.0 * obbi.extent[0] * obbi.e[0];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner += 2.0 * obbi.extent[1] * obbi.e[1];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
        corner -= 2.0 * obbi.extent[0] * obbi.e[0];
        fprintf(sortie_vtk, "%e %e %e\n", corner.x, corner.y, corner.z);
    }

    fprintf(sortie_vtk, "LINES %li %li\n", 6 * box.Particles.size(), 24 * box.Particles.size());

    long indexParticle = 0;

    for (size_t i = 0; i < box.Particles.size(); i++)
    {
        fprintf(sortie_vtk, "5 %li %li  %li %li %li\n", indexParticle, indexParticle + 1, indexParticle + 2, indexParticle + 3, indexParticle);
        fprintf(sortie_vtk, "5 %li %li  %li %li %li\n", indexParticle + 4, indexParticle + 5, indexParticle + 6, indexParticle + 7, indexParticle + 4);

        fprintf(sortie_vtk, "2 %li %li\n", indexParticle, indexParticle + 4);
        fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 1, indexParticle + 5);
        fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 2, indexParticle + 6);
        fprintf(sortie_vtk, "2 %li %li\n", indexParticle + 3, indexParticle + 7);
        indexParticle += 8;
    }

    fprintf(sortie_vtk, "CELL_DATA %li \n", 6 * box.Particles.size());
    fprintf(sortie_vtk, "SCALARS Id int 1\n");
    fprintf(sortie_vtk, "LOOKUP_TABLE default\n");
    for (size_t i = 0; i < box.Particles.size(); i++)
    {
        for (int j = 0; j < 6; j++)
        {
            fprintf(sortie_vtk, "%li\n", i);
        }
    }

    fclose(sortie_vtk);
}


void nbConfToVtk()
{
    std::string nomfich;
    char nom[25];
    struct dirent *lecture;
    DIR *rep;
    rep = opendir(".");

    int number;

    tabConf.clear();

    while ((lecture = readdir(rep)))
    {
        strcpy(nom, lecture->d_name);

        nomfich = nom;

        if (nomfich[0] == 'c' && nomfich[1] == 'o' && nomfich[2] == 'n' && nomfich[3] == 'f')
        {
            std::string temp;

            for (unsigned int i = 0; i < nomfich.size(); i++)
            {
                if (isdigit(nomfich[i]))
                {
                    temp += nomfich[i];
                }
            }
            istringstream stream(temp);
            stream >> number;
            tabConf.push_back(number);
        }
    }
    closedir(rep);

    // tirer :
    for (size_t i = 0; i < tabConf.size(); i++)
    {
        for (size_t j = i + 1; j < tabConf.size(); j++)
        {
            if (tabConf[j] < tabConf[i])
            {
                int aux = tabConf[i];
                tabConf[i] = tabConf[j];
                tabConf[j] = aux;
            }
        }
    }
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char *argv[])
{
    box.setInteractive(true);

    size_t ifirst = 0;
    size_t ilast = 0;
    size_t istep = 0;

    try
    {
        TCLAP::CmdLine cmd("VTK maker of Rockable simulations", ' ', "0.3");

        ifirst = std::atoi(argv[1]);
        ilast = std::atoi(argv[2]);
        istep = std::atoi(argv[3]);
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
    }

    // nbConfToVtk();

    for (size_t i = ifirst; i <= ilast; i += istep)
    {
        box.clearMemory();
        tryToReadConf(i);
        box.computeAABB();
        box.System.read();
        confNum = box.iconf;

        if (box.Particles.empty())
        {
            std::cerr << "No particles!" << std::endl;
        }
        else
        {
            SeparateParticlesByType();
            if (Sphers.size() > 0)
                writeVTKSPHER(confNum);
            if (Polyrs.size() > 0)
                writeVTKPOLYR(confNum);
            if (Joncs.size() > 0)
                writeVTKJONCx(confNum);
            // writeVTKContactsShperes(confNum);
            writeVTKContactsPolys(confNum);
            writeVTKOBB(confNum);
        }
    }
    return 0;
}