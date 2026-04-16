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

#include "factory.hpp"
#include "kwParser.hpp"

#include "PostProcessor_Network.hpp"
#include "Rockable.hpp"

#ifndef PERIODIC_BOX
#define PERIODIC_BOX
#endif

static Registrar<PostProcessor, Network> registrar("Network");

Network::Network() {}

bool Network::iscolinearity(const vec3r &v0, const vec3r &v1, const vec3r &v2)
{
  vec3r v01 = v1 - v0;
  v01.normalized();
  vec3r v02 = v2 - v0;
  v02.normalized();
  double cos0 = dot(v01, v02);
  double eps = 1.e-10;
  return ((1. - fabs(cos0)) < eps);
}

void Network::read(std::istream &is)
{
  kwParser parser;
  parser.parse(is);
}

void Network::init()
{
}

void Network::end()
{
}

void Network::exec()
{

  contacts.clear();
  box->System.read();

  box->activeInteractions.clear();

  // first construct the activeInteractions vector
  for (size_t k = 0; k < box->Interactions.size(); ++k)
  {
    for (auto it = box->Interactions[k].begin(); it != box->Interactions[k].end(); ++it)
    {
      Interaction *I = const_cast<Interaction *>(std::addressof(*it));
      if (it->dn < 0.0)
      {
        box->activeInteractions.push_back(I);
      }
    }
  }

  for (size_t k = 0; k < box->activeInteractions.size(); ++k)
  {
    Interaction *I = box->activeInteractions[k];
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

  std::sort(box->activeInteractions.begin(), box->activeInteractions.end(), std::less<Interaction *>());

  // now we construct the contacts vector
  size_t ncont = 0; // index related to contacts vector

  size_t i = box->activeInteractions[0]->i;
  size_t j = box->activeInteractions[0]->j;

  size_t ibefore = i;
  size_t jbefore = j;

  contacts.push_back(new Contact);
  contacts[ncont]->i = i;
  contacts[ncont]->j = j;

  for (size_t k = 0; k < box->activeInteractions.size(); ++k)
  {
    Interaction *I = box->activeInteractions[k];

    i = I->i;
    j = I->j;

    if (i != ibefore || j != jbefore) // i or j changes, so we add new Contact to the contacts vector
    {
      contacts.push_back(new Contact);
      ncont++; // we increment ncont after push_back
      contacts[ncont]->i = i;
      contacts[ncont]->j = j;
      contacts[ncont]->l = I->lji;            // the branch vector form i to j
      contacts[ncont]->subinter.push_back(I); // we add the interaction to the subinter
      ibefore = i;
      jbefore = j;
    }
    else // i and j doesn't change
    {
      contacts[ncont]->subinter.push_back(I); // here i and j dont changes so we just add the interaction to the subinter
    }
  }

  // now we calculate the fn and ft vector (vectorial sum)
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

  }

  // now the type of contact

  size_t sc = 0, dc = 0, tc = 0;

  for (size_t c = 0; c < contacts.size(); ++c)
  {

    Contact *cij = contacts[c];
    size_t nbc = cij->subinter.size(); // number of 'rockable' contacts between i and j

    if (nbc > 4)
    {
      tc++;
      cij->type = t_c;
    } // if number of contacts is bigger than 4 it is a face face contact
    else // we have to check to define a contact type
    {
      vec3r pos0, pos1, pos2, pos3, pos;
      double eps;

      switch (nbc)
      {
      case 0:
        std::cout << "No Contact Type !!" << std::endl;
        break;
      case 1:
        sc++;
        cij->type = s_c;
        break;
      case 2:
        pos0 = cij->subinter[0]->posI;
        pos1 = cij->subinter[1]->posI;
        eps = 1.e-04 * box->Particles[0].MinskowskiRadius();
        pos = pos1 - pos0;
        pos.round(box->periodicity);
        if (norm(pos) < eps)
        {
          sc++;
          cij->type = s_c;
        }
        else
        {
          dc++;
          cij->type = d_c;
        }
        break;
      case 3:
        pos0 = cij->subinter[0]->posI;
        pos1 = cij->subinter[1]->posI;
        pos2 = cij->subinter[2]->posI;
        if (iscolinearity(pos0, pos1, pos2))
        {
          dc++;
          cij->type = d_c;
        }
        else
        {
          tc++;
          cij->type = t_c;
        }
        break;
      case 4:
        pos0 = cij->subinter[0]->posI;
        pos1 = cij->subinter[1]->posI;
        pos2 = cij->subinter[2]->posI;
        pos3 = cij->subinter[3]->posI;
        if (iscolinearity(pos0, pos1, pos2) && iscolinearity(pos1, pos2, pos3))
        {
          dc++;
          cij->type = d_c;
        }
        else
        {
          tc++;
          cij->type = t_c;
        }
        break;

      default:
        std::cout << "Contact Type Probleme" << std::endl;
        break;
      }
    }
  }



}
