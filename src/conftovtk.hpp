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

#ifndef CONFTOVTK_HPP
#define CONFTOVTK_HPP

#include "CmdLine.h"
#include "Rockable.hpp"
#include "message.hpp"

#include <string.h>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>

#define PLANx "PLANx"
#define SPHER "SPHER"
#define POLYR "POLYR"

#define s_c 0 // to define contact type
#define d_c 1
#define t_c 2

//@dv
struct Contact
{
  size_t i;   // index particle i
  size_t j;   // index particle j
  vec3r fn;   // normal force
  vec3r ft;   // tangential force
  vec3r l;    // branch vector
  vec3r n;    // normal unit vector
  vec3r t;    // tangential unit vector
  vec3r posI; // position of total contact
  vec3r posJ;
  size_t type;                         // contact type
  std::vector<Interaction *> subinter; // interactions betwween i and j;
};

std::vector<Contact *> contacts; // all contacts vector

Rockable box;
std::vector<Particle *> Sphers;
std::vector<Particle *> Joncs;
std::vector<Particle *> Polyrs;

std::vector<size_t> SphersId;
std::vector<size_t> JoncsId;
std::vector<size_t> PolyrsId;

int confNum = 0;
std::vector<int> tabConf;
size_t complexityNumber = 0; // it says how the sample will be long to display
bool tryToReadConf(int num);
void SeparateParticlesByType();
void nbConfToVtk();           // nombre de configurations
void writeVTKSPHER(int num);  // write spheres
void writeVTKPOLYR(int num);  // write polyhedra
void writeVTKJONCx(int num);  // write elongated
void writeVTKOBB(int num);    // write Obbs box

std::string TruncateLongstring(std::string const str, unsigned int maxLength);

bool iscolinearity(const vec3r &v0, const vec3r &v1, const vec3r &v2);

#endif /* end of include guard: SEE_HPP_E29BD15E */
