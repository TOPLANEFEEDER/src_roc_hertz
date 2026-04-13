//
//  tessToshp.hpp
//  rockable
//
//  Created by amarsid on 25/09/2019.
//  Copyright © 2019 amarsid. All rights reserved.
//

#ifndef tessToshp_hpp
#define tessToshp_hpp

#define CONF_VERSION_DATE "29-11-2018"

#include "Rockable.hpp"
#include "polyhToolCut.hpp"

template <class myType>
int getIndex ( std::vector<myType> v, myType b) {
 	auto it = find(v.begin(), v.end(), b);
 	if (it != v.end()) return int(it - v.begin());
 	return -1;
}

struct plan_
{
	vec3r pos;
	vec3r normal;
	double d;
};

struct vertex_
{
	int id;
	vec3r v;
	bool operator ==(const vertex_& a){ return std::abs(id) == std::abs(a.id);}
	bool operator <(const vertex_& a){ return std::abs(id) < std::abs(a.id);}
};

struct edge_
{
	int id;
	int ver_1;
	int ver_2;
	bool operator ==(const edge_& a){ return std::abs(id) == std::abs(a.id);}
	bool operator <(const edge_& a){ return std::abs(id) < std::abs(a.id);}
};

struct face_
{
	int id;
	std::vector<int> ver_;
	std::vector<int> edge_;
	double face_eq_d, face_eq_a, face_eq_b, face_eq_c;
	bool operator ==(const face_& a){ return std::abs(id) == std::abs(a.id);}
	bool operator <(const face_& a){ return std::abs(id) < std::abs(a.id);}
};

struct polyhedron_
{
	int id;
	std::vector<int> face_;
	bool operator ==(const polyhedron_& a){ return std::abs(id) == std::abs(a.id);}
};

struct shape_
{
	std::vector<vec3r> v ;
    std::vector<std::pair<size_t,size_t>> e;
    std::vector<std::vector<size_t>> f;
};

class tessToshp : public Rockable
{

public:
	tessToshp(){};
	virtual ~tessToshp(){};
	void createShpere();
	void addPlan(double theta,double psi);
	void createTess();
	void readTess();
	void extractShapes();
	void validate_shape_faces(size_t num_shp);
	void massProperties();
	void writeShapes();
	void loadConfNeper(const char* name);
	void saveConfNeper();
	bool all_vertexes_in_shape(Shape & shape, int & iv_out);
	vec3r mass_center(Shape & shp);
	void exec();
	void create_shapes();
	void CalculateAreas(Shape & shp);
	
private:
	std::vector<vec3r> position_;
	std::vector<vertex_> vertex;
	std::vector<edge_> edge;
	std::vector<face_> face;
	std::vector<polyhedron_> polyhedron;
	std::vector<plan_> plans;
	std::vector<std::vector<plan_>> face_plans;
	std::vector<std::vector<plan>> shape_plans;
	double R; //rayon de la shepere;
	int nbP; //nombre de plans
	double radius;
	// les noms de variables sont définis dans la documentation de Neper
	int number_of_cells;
	double min_norm_edge;
	double max_norm_edge;
};

#endif /* tessToshp_hpp */
