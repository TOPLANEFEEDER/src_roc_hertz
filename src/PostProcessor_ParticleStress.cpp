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

#include "PostProcessor_ParticleStress.hpp"
#include "Rockable.hpp"

static Registrar<PostProcessor, ParticleStress> registrar("ParticleStress");

ParticleStress::ParticleStress() {}

void ParticleStress::read(std::istream &is)
{
  kwParser parser;
  parser.kwMap["Volume"] = __GET__(is, Volume);
  parser.kwMap["ConfVolumes"] = __DO__(is)
  {
    size_t nb;
    is >> nb;
    for (size_t i = 0; i < nb; i++)
    {
      int iconf;
      double v;
      is >> iconf >> v;
      ConfVolumes[iconf] = v;
    }
  };
  parser.parse(is);
}

void ParticleStress::setVolume()
{
  box->System.read();
  Volume = box->periodicity.x * box->periodicity.y * box->periodicity.z;
}

void ParticleStress::init()
{
}

void ParticleStress::end()
{
}

void ParticleStress::exec()
{

  setVolume();

  mat9r Mtotal;
  for (size_t ibody = 0; ibody < box->Particles.size(); ibody++)
  {
    for (auto it = box->Interactions[ibody].begin(); it != box->Interactions[ibody].end(); ++it)
    {

      Interaction *I = const_cast<Interaction *>(std::addressof(*it));

      vec3r rji = I->lji;
      vec3r f = (I->fn * I->n) + I->ft; // the sign - is for obtaining positive M for compression

      if (I->dn <= 0.)
      {
        double sxx = f.x * rji.x;
        double sxy = f.x * rji.y;
        double sxz = f.x * rji.z;
        double syx = f.y * rji.x;
        double syy = f.y * rji.y;
        double syz = f.y * rji.z;
        double szx = f.z * rji.x;
        double szy = f.z * rji.y;
        double szz = f.z * rji.z;

        Mtotal.xx += sxx;
        Mtotal.xy += sxy;
        Mtotal.xz += sxz;
        Mtotal.yx += syx;
        Mtotal.yy += syy;
        Mtotal.yz += syz;
        Mtotal.zx += szx;
        Mtotal.zy += szy;
        Mtotal.zz += szz;
      }
    }
  }

  /*for (size_t ibody = 0; ibody < box->Particles.size(); ibody++)
  {
    double m = box->Particles[ibody].mass;
    vec3r v = box->Particles[ibody].vel;

    Mtotal.xx += m * v.x * v.x;
    Mtotal.xy += m * v.x * v.y;
    Mtotal.xz += m * v.x * v.z;
    Mtotal.yx += m * v.y * v.x;
    Mtotal.yy += m * v.y * v.y;
    Mtotal.yz += m * v.y * v.z;
    Mtotal.zx += m * v.z * v.x;
    Mtotal.zy += m * v.z * v.y;
    Mtotal.zz += m * v.z * v.z;
  }*/

  double Vtot = 1.0;

  Vtot = Volume;
  __SHOW(box->iconf);
  __SHOW(Vtot);

  std::cout << "total: " << Mtotal / Vtot << std::endl;
  std::ofstream out;
  out.open("TotalStress.txt", std::ofstream::out | std::ofstream::app);
  out << std::scientific << std::setprecision(6);
  out << box->t << ' ' << Mtotal / Vtot << std::endl;
  out.close();
}
