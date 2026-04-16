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

#ifndef POSTPROCESSOR_NETWORK_CPP
#define POSTPROCESSOR_NETWORK_CPP

#define s_c 0 // to define contact type
#define d_c 1
#define t_c 2

#include "PostProcessor.hpp"
#include "Interaction.hpp"
#include "vec3.hpp"

struct Contact
{
  size_t i;                            // index particle i
  size_t j;                            // index particle j
  vec3r fn;                            // normal force
  vec3r ft;                            // tangential force
  vec3r l;                             // branch vector
  vec3r n;                             // normal unit vector
  vec3r t;                             // tangential unit vector
  size_t type;                         // contact type
  std::vector<Interaction *> subinter; // interactions betwween i and j;
};

class Network : public PostProcessor
{
public:
  Network();
  void init();
  bool iscolinearity(const vec3r &v0, const vec3r &v1, const vec3r &v2);
  void end();
  void read(std::istream &is);
  void exec();

private:
  std::vector<Contact *> contacts; // all contacts vector
};

#endif /* end of include guard: POSTPROCESSOR_NETWORK_CPP */