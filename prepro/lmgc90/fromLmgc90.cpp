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

#include "CmdLine.h"
#include "lmgc90.hpp"

int main(int argc, char const* argv[]) {
  std::string confFileName;
  int nbThreads = 1;

  try {
    TCLAP::CmdLine cmd("This is the main runner application for Rockable", ' ', "0.3");
    TCLAP::UnlabeledValueArg<std::string> nameArg("inputlmgc90", "Name of the conf-file", true, "conf0", "conf-file");
    TCLAP::ValueArg<int> nbThreadsArg("j", "nbThreads", "Number of threads to be used", false, 1, "int");

    cmd.add(nameArg);
    cmd.add(nbThreadsArg);

    cmd.parse(argc, argv);

    confFileName = "inputlmgc90.txt";//nameArg.getValue();
    nbThreads = nbThreadsArg.getValue();
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for argument " << e.argId() << std::endl;
  }

  std::cout << "No acceleration" << std::endl;
  if (nbThreads != 1) std::cout << "nbThreads = 1\n";

  lmgc90  box ;
  box.showBanner();
  box.start();
  box.loadConfLMGC90(confFileName.c_str());  
  box.readFromLMGC90();
  box.defineALL();
  if(box.isdrumDefined())
  {
      box.InitDrum();
      if(box.drum_type()==CYL_T0) box.defineDrum();
      if(box.drum_type()==CYL_T1) box.defineDrums();
      if(box.drum_type()==CYL_T2) box.defineCylinder();
  }
  box.saveConfLMGC90();
  box.writeLmgcShape();
  std::cout << msg::info() << "LMGC TO ROCKABLE NORMALLY STOPPED" << msg::normal() << std::endl;
  return 0;
}
