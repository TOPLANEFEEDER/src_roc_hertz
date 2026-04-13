//
//  lmgc90.hpp
//  rockable
//
//  Created by amarsid on 25/09/2019.
//  Copyright © 2019 amarsid. All rights reserved.
//

#ifndef lmgc90_hpp
#define lmgc90_hpp

#define CONF_VERSION_DATE "29-11-2018"

#include "Cylinder.hpp"
#include "Rockable.hpp"

#define PLANx "PLANx"
#define SPHER "SPHER"
#define POLYR "POLYR"

#define CYL_T0 0
#define CYL_T1 1
#define CYL_T2 2

struct shapeProp {
    size_t ishape;
    double avrd;
    double I1;
    double I2;
    double I3;
    double density;
    double volume;
};

class lmgc90 : public Rockable {

public:
    lmgc90(){};
    virtual ~lmgc90(){};

    void readBULK_BEHAV();
    void readBODIES();
    void readDOF();
    void readTACT_BEHAV();
    void defineALL();

    void readFromLMGC90();
    void loadConfLMGC90(const char *name); ///< Load a configuration to use lmgc90
    void saveConfLMGC90();                 ///< save a configuration to use lmgc90

    void readPLANx(std::ifstream &in, Particle &P, std::string &color, std::string &material, double &avrd, double &I1, double &I2, double &I3);
    void readSPHER(std::ifstream &in, Particle &P, std::string &color);
    void readPOLYR(std::ifstream &in, Particle &P, std::string &color, std::string &material, double &avrd, double &I1, double &I2, double &I3);

    Shape *createSpherShape(double radius);
    Shape *createPlanxShape(double &axe1, double &axe2, double &axe3);
    void addPolyrShape(Shape *shape);
    void writeLmgcShape();

    void setVolume_Inirtia_OBB(Shape *shape, shapeProp &shp_pro);

    bool existSpherShape; //pour les shperes
    int existPlanxShape(double &demiaxe1, double &demiaxe2, double &demiaxe3);
    int existPolyrShape(Shape *shape, shapeProp &shapepro);

    size_t shape_by_name(std::string &name);

    std::string truncateLongString(std::string &str, unsigned int maxLength);
    std::string convertToString(char *a, size_t size);

    void exacteSize(Particle &P);
    void exacteSize(Shape *shape);
    void updateshapeId();
    void start();

    int netoyagePolyr(Shape *shape);
    void reduceSizeShape(Shape *shape, double &radius);
    bool siTowFacesConnected(std::vector<size_t> face1, std::vector<size_t> face2);
    std::pair<size_t, size_t> commonEdgeFacesConnected(std::vector<size_t> face1, std::vector<size_t> face2);

    unsigned int extractnumber(std::string &str);
    void removeEgde(std::vector<std::pair<size_t, size_t>> &vedge, std::pair<size_t, size_t> edge);
    void addEgde(std::vector<std::pair<size_t, size_t>> &vedge, std::pair<size_t, size_t> edge);
    void concatener_deux_faces(std::vector<std::vector<size_t>> &face, size_t fmax, size_t fmin, std::pair<size_t, size_t> edge);
    vec3r get_normal_face(size_t &f, Shape *shape);
    void removeFace(std::vector<std::vector<size_t>> &face, size_t f1);

    void writeCorresondingSpheres(std::ofstream &out);

    void defineDrum();
    void defineDrums();
    void defineCylinder();
    void InitDrum();
    bool isdrumDefined() { return isdrum; };
    size_t drum_type() { return type_drum; };
    void addDrum();
    bool isInCylinder(Particle *P);
    bool isInCylinder(vec3r &v);
    bool isInCylinder(vec3r &v, double &dmax);

private:
    int iShapes_;
    double rInput;
    vec3r aabbextent;

    //pour extraire les lois, les materiaux
    std::map<std::pair<std::string, std::string>, size_t> groupes_;
    std::map<std::string, double> law_fric_; //pour chaque loi on definie une friction
    std::vector<std::string> isqc_tab_;
    std::vector<std::pair<std::string, std::string>> cans_type_color_;
    std::vector<std::pair<std::string, std::string>> ants_type_color_;

    std::vector<std::string> colors_;         //colors definie les groupes
    std::vector<std::string> materials_;      //material definie les densites
    std::map<std::string, double> densities_; //pour chaque materiaux on a une densite

    std::map<Particle *, double> Hmax_; //homoticy initiale par particule;
    std::vector<shapeProp> shape_pro;   //garder les proprietes de shape;

    std::vector<Particle *> insideCyl; //Particles inside Cylinder

    //Les identifiants POLYR
    std::vector<vec3r> id_shape_PLANx_;
    size_t id_shape_POLYR_;
    std::vector<double> avrd_polyr_;
    std::vector<std::string> index_shapes_;

    //for drum
    bool isdrum;
    size_t type_drum, group_drum;
    bool correspondingSphere;
    Shape drumShp;
    std::vector<Particle> drums;
    double fric_drum;
    size_t resolution;
    double radius_drum;
    double kn, kt, kr, en, mur;
    Cylinder drum_cyl;
};

#endif /* lmgc90_hpp */
