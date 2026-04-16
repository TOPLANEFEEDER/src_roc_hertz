//
//  Lmgc90.cpp
//  rockable
//
//  Created by amarsid on 25/09/2019.
//  Copyright © 2019 amarsid. All rights reserved.
//

#include "Lmgc90.hpp"
#include "polyhTool.hpp"
#include <cstdlib>

//IMPORTANT :  orientation.set_rot_matrix(VP.c_mtx()); pour orienter à partir du lmgc 90

std::string Lmgc90::truncateLongString(std::string & str, unsigned int maxLength)
{
    if(str.length() > maxLength)
        return  str.substr(0,maxLength);
    return str;
}

std::string Lmgc90::convertToString(char* a, size_t size) 
{ 
    size_t i; 
    std::string s = ""; 
    for (i = 0; i < size; i++) { 
        s = s + a[i]; 
    } 
    return s; 
} 

void Lmgc90::start()
{
    //des variables à initialiser
    iShapes_ =0;  
    rInput = 0.;
    aabbextent = vec3r();
    index_shapes_.clear();
    id_shape_PLANx_.clear();
    id_shape_POLYR_=0;
    avrd_polyr_.clear();

    //une sphere n a qu un seul shape :-)
    existSpherShape = false;

    //les vectors de la classe rockable
    Shapes.clear();
    shapeId.clear();
    Particles.clear();

    //parametres contact
    kn = kt = kr = mur = en =0.;
}

void Lmgc90::readBULK_BEHAV()
{
    std::string fin_ = "./DATBOX/BULK_BEHAV.DAT";
    std::ifstream in(fin_.c_str());
    std::string line_,stemp_;
    
    if(in)
    {
        vec3r & g = gravity;

        do std::getline(in, line_); while (truncateLongString(line_,6)!=std::string("$gravy"));// (line!=string("$gravy  "));
        std::getline(in, line_);
        sscanf(line_.c_str(),"                   grv1=%le  grv2=%le  grv3=%le\n",&g.x,&g.y,&g.z);
        
        while(!in.eof())
        {
            std::string material_;
            double density_;
            
            std::getline(in, line_); line_ = truncateLongString(line_,6);
            if(line_==std::string("$behav"))
            {
                in>>material_>>stemp_>>stemp_>>density_;
                densities_[material_] = density_;
            }
        }
        in.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fin_.c_str() <<  std::endl;
    }
}

void Lmgc90::readTACT_BEHAV()
{
    dVerlet = 1e+30;
    DVerlet = -1e+30;

    std::string fin_ = "./DATBOX/TACT_BEHAV.DAT";
    std::ifstream in(fin_.c_str());
    std::string line_,stemp_,type0_,type1_,iqsc_,color0_,color1_;
    double friction_=0.;
    double alert_=0.;
    if(in)
    {   
        while(!in.eof())
        {
            std::getline(in, line_); line_ = truncateLongString(line_,6);

            if(line_==std::string("$behav"))
            {
                in>>iqsc_>>stemp_>>stemp_>>friction_;
                law_fric_[iqsc_] = friction_;
            }

            if(line_==std::string("$seety"))
            {
                std::getline(in, line_);
                in>>stemp_>>type0_>>color0_>>iqsc_>>stemp_>>type1_>>color1_>>alert_;
                
                dVerlet = std::min(dVerlet,alert_);
                DVerlet = std::max(DVerlet,alert_);

                std::pair<std::string,std::string> can_type_color_;
                std::pair<std::string,std::string> ant_type_color_;
                
                can_type_color_.first = type0_;
                can_type_color_.second = color0_;
                ant_type_color_.first = type1_;
                ant_type_color_.second = color1_;
                cans_type_color_.push_back(can_type_color_);
                ants_type_color_.push_back(ant_type_color_);
                isqc_tab_.push_back(iqsc_);
            }
        }
        in.close();
        /*
        //verification
        std::cout<<"------------law_fric_-----------------"<<std::endl;
        std::map<std::string,double>::iterator it;
        for(it = law_fric_.begin(); it!=law_fric_.end(); ++it)
        {
            std::cout << it->first <<" "<<it->second << std::endl; 
        }
        std::cout<<"-----------------------------------------"<<std::endl;
            
        //verification
        std::cout<<"------------cans_type_color_-----------------"<<std::endl;
        for(size_t i = 0; i<cans_type_color_.size(); ++i)
        {
            std::cout << cans_type_color_[i].first <<" "<< cans_type_color_[i].second << std::endl; 
        }
        std::cout<<"-----------------------------------------"<<std::endl;

        std::cout<<"------------ants_type_color_-----------------"<<std::endl;
        for(size_t i = 0; i<ants_type_color_.size(); ++i)
        {
            std::cout << ants_type_color_[i].first <<" "<< ants_type_color_[i].second << std::endl; 
        }
        std::cout<<"-----------------------------------------"<<std::endl;
        
        std::cout<<"------------isqc_tab_-----------------"<<std::endl;
        for(size_t i = 0; i<isqc_tab_.size(); ++i)
        {
            std::cout << isqc_tab_[i] << std::endl; 
        }
        std::cout<<"-----------------------------------------"<<std::endl;
        */
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fin_.c_str() <<  std::endl;
    }
}

void Lmgc90::readBODIES()
{
    //on lit le fichier qui contient les particules
    std::string line_="";


    std::string type_="";
    Particle P_ = Particle();
    double avrd_,I1_,I2_,I3_;
    size_t i_groupe_ = 0;

    std::string fin_ = "./DATBOX/BODIES.DAT";
    std::ifstream in(fin_.c_str());
    
    if(in)
    {
        while(!in.eof())
        {
            std::getline(in, line_); line_ = truncateLongString(line_,6);
            if(line_==std::string("$bdyty"))
            {
                std::string  strTemp="";
                std::string  material_="";
                std::string  blmty = "$blmty";
                std::string  nodty = "$nodty";
                std::string  tacty = "$tacty";

                int itemp_=0;

                do std::getline(in, line_); while(line_.find(blmty) == std::string::npos);//(truncateLongString(line_,6)!=blmty);
                in>>strTemp>>itemp_>>strTemp>>material_>>strTemp>>avrd_;
                std::getline(in, line_);
                std::getline(in, line_);
                std::sscanf(line_.c_str()," I1  = %le  I2  = %le  I3  = %le\n",&I1_,&I2_,&I3_);
                do std::getline(in, line_); while(line_.find(nodty) == std::string::npos);//while (truncateLongString(line_,6)!=nodty);
                
                in>>strTemp>>itemp_;
                std::getline(in, line_);
                std::sscanf(line_.c_str(),"                coo1=%le  coo2=%le  coo3=%le\n",&P_.pos.x,&P_.pos.y,&P_.pos.z);
                do std::getline(in, line_); while(line_.find(tacty) == std::string::npos);//while (truncateLongString(line_,6)!=tacty);
                
                in>>type_;

                std::string color;

                if(type_==SPHER)readSPHER(in,P_,color);
                if(type_==POLYR)readPOLYR(in,P_,color,material_,avrd_,I1_,I2_,I3_);
                if(type_==PLANx)readPLANx(in,P_,color,material_,avrd_,I1_,I2_,I3_);

                materials_.push_back(material_);
                colors_.push_back(color);

                std::pair<std::string,std::string> pair_;
                pair_.first  = type_;
                pair_.second = color;

                if(groupes_.insert(make_pair(pair_,i_groupe_)).second)
                {
                    //std::cout<<"pair dans readBodies "<<pair_.first<<" "<<pair_.second<<std::endl;
                    //printf("material   %s    density    %e\n",material_,densities_[material_]);
                    //std::cout<<"material   "<<material_<<"    density    "<<densities_[material_]<<std::endl;
                    groupes_[pair_] = i_groupe_;
                    properties.set(idDensity, i_groupe_, densities_[material_]);
                    P_.group = i_groupe_;
                    i_groupe_++;
                }
                else
                {P_.group = groupes_[pair_];}

                Particles.push_back(P_);
           }
        }
        in.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fin_.c_str() <<  std::endl;
    }
}

void Lmgc90::readDOF()
{
    std::string line_;
    
    //ce fichier contient les rotations des particules 
    std::string fin_ = "./DATBOX/DOF.INI";
    std::ifstream in(fin_.c_str());
    if(in)
    {
        unsigned int index_p_  = 0;
        double M0_[9];
        double M1_[9];
        double M2_[9];
        mat9r m1_m0_,m2_m0_;
        while(!in.eof())
        {
            std::getline(in, line_); line_ = truncateLongString(line_,6);
            if(line_==std::string("$bdyty"))
            {
                std::string type_;
                Particle & P = Particles[index_p_];

                do std::getline(in, line_); while (truncateLongString(line_,6)!=std::string("$nodty"));
                for(int i=0;i<5;i++) std::getline(in, line_);
                std::sscanf(line_.c_str(),"                             a(1)=%le  a(2)=%le  a(3)=%le \n",&M0_[0],&M0_[3],&M0_[6]);
                std::getline(in, line_);
                std::sscanf(line_.c_str(),"                             b(1)=%le  b(2)=%le  b(3)=%le \n",&M0_[1],&M0_[4],&M0_[7]);
                std::getline(in, line_);
                std::sscanf(line_.c_str(),"                             c(1)=%le  c(2)=%le  c(3)=%le \n",&M0_[2],&M0_[5],&M0_[8]);

                //Choix du bon quaternion car la fonction set_rot_matrix peut donner soit +Q ou -Q (congugate(Q))
                quat Q,Qi;
                Q.set_rot_matrix(M0_);  

                Qi = Q;
                Qi.conjugate();

                Q.get_rot_matrix(M1_);
                Qi.get_rot_matrix(M2_);

                double fabs_m1_m0_ = 0.;
                double fabs_m2_m0_ = 0.;

                for(size_t i=0;i<9;i++)
                {
                    fabs_m1_m0_ += fabs(M1_[i] - M0_[i]);
                    fabs_m2_m0_ += fabs(M2_[i] - M0_[i]);
                }

                if(fabs_m1_m0_<fabs_m2_m0_) P.Q = Q;
                else P.Q = Qi;

                index_p_++;

                /*
                //test Compatibilite
                std::cout<<"M0"<<std::endl;
                std::cout<<M_[0]<<" "<<M_[1]<<" "<<M_[2]<<std::endl;
                std::cout<<M_[3]<<" "<<M_[4]<<" "<<M_[5]<<std::endl;
                std::cout<<M_[6]<<" "<<M_[7]<<" "<<M_[8]<<std::endl;
                P.Q.get_rot_matrix(M_);
                std::cout<<"M1"<<std::endl;
                std::cout<<M_[0]<<" "<<M_[1]<<" "<<M_[2]<<std::endl;
                std::cout<<M_[3]<<" "<<M_[4]<<" "<<M_[5]<<std::endl;
                std::cout<<M_[6]<<" "<<M_[7]<<" "<<M_[8]<<std::endl;
                getchar();
                */
            }
        }
        in.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fin_.c_str() <<  std::endl;
    }
}

void Lmgc90::writeLmgcShape()
{
    std::string fout_ = "./shape.shp";
    std::ofstream out(fout_.c_str());
    
    shapeFile = fout_.c_str();

    if(out)
    {
        for(size_t s = 0; s < Shapes.size(); ++s) Shapes[s].write(out);
        if(correspondingSphere)
        {
            //writeCorresondingSpheres(out);
        }
        out.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fout_.c_str() <<  std::endl;
    }
}

void Lmgc90::writeCorresondingSpheres(std::ofstream &out)
{
    for(size_t i=0;i<Shapes.size();i++)
    if(Shapes[i].vertex.size()>1)
    {
        Shape shp = Shapes[i];

        double max_extent = 0.;
        max_extent = std::max(shp.obb.extent.x,max_extent);
        max_extent = std::max(shp.obb.extent.y,max_extent);
        max_extent = std::max(shp.obb.extent.z,max_extent);
        
        Shape s = createSpherShape(max_extent);
        s.name = "Sphere" + shp.name;
        s.write(out);
    }
}

void Lmgc90::defineALL()
{
    size_t g1;
    size_t g2;

/*
    std::cout<<"-----------Verification des groupes------------"<<std::endl;
    std::map<std::pair<std::string,std::string>,size_t>::iterator it;
    for(it = groupes_.begin(); it!=groupes_.end(); ++it)
    {
        std::cout << it->first.first <<" "<<it->first.second<<" "<<it->second << std::endl; 
    }
    std::cout<<"-----------------------------------------"<<std::endl;

    std::cout<<"-----------Construction des groupes------------"<<std::endl;
    for(size_t i=0;i<cans_type_color_.size();i++)
    {
        g1 = groupes_[cans_type_color_[i]];
        g2 = groupes_[ants_type_color_[i]];
        std::cout<<"g1 "<<g1<<" "<<g2<<std::endl;
    }
    std::cout<<"-----------------------------------------------"<<std::endl;
*/
    //definir les groupes
    for(size_t i=0;i<cans_type_color_.size();i++)
    {
        g1 = groupes_[cans_type_color_[i]];
        g2 = groupes_[ants_type_color_[i]];
        dataTable.set (idMuContact, g1, g2, law_fric_[isqc_tab_[i]]);
        dataTable.set (idKnContact, g1, g2, kn);
        dataTable.set (idKtContact, g1, g2, kt);
        dataTable.set (idKrContact, g1, g2, kr);
        dataTable.set (idMurContact, g1, g2, mur);
        dataTable.set (idEn2Contact, g1, g2, en);
    }

    /*
    //les densites
    std::cout<<"-----------Les densites------------"<<std::endl;
    for (size_t grp = 0; grp < properties.ngroup; grp++) {
        double density = properties.get(idDensity, grp);
            std::cout << "density " << grp << " " << density << std::endl;
    }
    std::cout<<"-----------------------------------------------"<<std::endl;
    */
}

void Lmgc90::readFromLMGC90()
{
    readBULK_BEHAV();
    readTACT_BEHAV();    
    readBODIES();
    readDOF();
    /*
    std::cout<<"verification des shapes"<<std::endl;
    std::cout<<"Nombre shape "<<Shapes.size()<<std::endl;
    for(size_t i = 0;i<Shapes.size();i++)
    {
        std::cout<<"nom shape "<<Shapes[i].name<<" Adress "<<&(Shapes[i])<<std::endl;
    }
    std::cout<<"verification des index_shapes_ size Particles"<<Particles.size()<<std::endl;
    for(size_t i = 1725;i<Particles.size();i++)
    {
            std::cout<<"nom shape "<<index_shapes_[i]<<" Adress "<<Particles[i].shape<<std::endl;
    }*/
    //definir les proprietes des shapes
    if(rInput==0.)
    {   
        rInput = 1e+30;
        for(size_t i=0;i<shape_pro.size();i++){rInput = std::min(rInput,shape_pro[i].avrd);}
        for(size_t i=0;i<shape_pro.size();i++){Shapes[i].radius = rInput/6.;}
    }
    else
    {
        for(size_t i=0;i<shape_pro.size();i++)
        Shapes[i].radius = rInput;
    }
    
    //je dis que le rayon des shape est egale au volume min / 10 ;
    for(size_t i=0;i<Shapes.size();i++)
    {   
        setVolume_Inirtia_OBB(Shapes[i],shape_pro[i]);    
    }

    updateshapeId();
    /*
    std::cout<<"verification de shapeId, taille "<<shapeId.size()<<std::endl;
    std::map<std::string,size_t>::iterator pos;
    for ( pos = shapeId.begin(); pos != shapeId.end(); ++pos)
    {
        std::cout<<"shapeId.frst "<<pos->first<<" shapeId.snd "<<pos->second<<std::endl;
    }
    */
    //une fois qu on a lu les shapes on les affecte aux particules 
    for(size_t i = 0;i<Particles.size();i++)
    {
        Particle & P = Particles[i];
        Particles[i].shape = &(Shapes[shapeId[index_shapes_[i]]]);
        Hmax_[&P] = P.homothety;
    }

    
    //maintenant on va mettre la taille exacte des particules
    /*
    for(size_t i = 0;i<Particles.size();i++)
    {
        exacteSize(Particles[i]);
    }*/


    for(size_t i = 0;i<Shapes.size();i++)
    {
        exacteSize(Shapes[i]);
    }
    

    /*
    for(size_t s = 0;s<Shapes.size();s++)
    reduceSizeShape(Shapes[s],Shapes[s].radius);
    */

    //FromDEM3DFaces();
    //ToDEM3DFaces();
    
    //un swap pour mettre les plan qui sont fixes au tout debutu
    for(size_t i = 0;i<Particles.size();i++)
    {
        for(size_t j = i+1;j<Particles.size();j++)
        {
            std::string str_i = Particles[i].shape->name;
            std::string str_j = Particles[j].shape->name;
            
            if((truncateLongString(str_j,5)==PLANx) && (truncateLongString(str_i,5)!=PLANx) ) std::swap(Particles[i],Particles[j]);
        }
    }
    
    if (Interactions.size() != Particles.size()) Interactions.resize(Particles.size());
    if (Interfaces.size() != Particles.size())   Interfaces.resize(Particles.size());

    double maxRadius = -1.e+30;
    double minRadius =  1.e+30;

    dVerlet = 1.e+30;
    DVerlet = -1.e+30;

    for(size_t s = 0; s < Shapes.size(); ++s)
    {
        dVerlet = std::min(dVerlet,Shapes[s].radius);
        DVerlet = std::max(DVerlet,Shapes[s].radius);
    }
    
    for(size_t s = 0; s < avrd_polyr_.size(); ++s)
    {
        minRadius = std::min(minRadius,avrd_polyr_[s]);
        maxRadius = std::max(maxRadius,avrd_polyr_[s]);
    }
    
    std::cout<<"maxRadius "<<maxRadius<<" minRadius "<<minRadius<<std::endl;

    dVerlet = std::min(DVerlet,minRadius);
    DVerlet = (minRadius + maxRadius)/2.;
}

void Lmgc90::setVolume_Inirtia_OBB(Shape & shape, shapeProp & shp_pro)
{
    AABB box(shape.vertex);   
    box.enlarge(shape.radius);
    shape.defineVertexConnectivity();
    shape.preCompDone='y';
    shape.volume = (4./3.)*Mth::pi*shp_pro.avrd*shp_pro.avrd*shp_pro.avrd;
    //Attention sur LMGC90 le I1,I2,I3 = I/density;
    shape.inertia_mass = vec3r(shp_pro.I1,shp_pro.I2,shp_pro.I3)/shape.volume;
    shape.obb.extent = (box.max-box.min)/2.;
    shape.obb.center = (box.max + box.min)/2.;
    shape.obb.e[0] = vec3r(1.,0.,0.);
    shape.obb.e[1] = vec3r(0.,1.,0.);
    shape.obb.e[2] = vec3r(0.,0.,1.);
    /*
    std::cout<<".................setVolume_Inirtia_OBB......................."<<std::endl;
    std::cout<<"shape name "<<shape.name<<std::endl;
    std::cout<<"shape volume "<<shape.volume<<std::endl;
    std::cout<<"shape density "<<density<<" mass "<<m<<std::endl;
    std::cout<<"I1 "<<I1<<" I2 "<<I2<<" I3 "<<I3<<std::endl;
    std::cout<<"............................................................"<<std::endl;
    */
    shape.clean();
    shape.defineVertexConnectivity();
}

void Lmgc90::readSPHER(std::ifstream & in, Particle & P, std::string & color)
{
    if(!existSpherShape){
        Shape shape = createSpherShape(1.);
        Shapes.push_back(shape);
        index_shapes_.push_back(shape.name);
        existSpherShape = true;
        iShapes_++;
    }
    else
    {
        index_shapes_.push_back(SPHER);
    }
    
    //homothety
    std::string stemp;
    int itemp;
    in>>itemp>>stemp>>color>>stemp>>P.homothety;
}

void Lmgc90::readPLANx(std::ifstream & in, Particle & P, std::string & color, std::string & material, double & avrd, double & I1, double & I2, double &I3)
{
    std::string line;
    int shapeId_=-1 ;
    double demiaxe1,demiaxe2,demiaxe3;
    int itemp;
    char color_[5];
    
    std::getline(in, line);
    
    std::sscanf(line.c_str(),"%d  color  %s  axe1=%le  axe2=%le  axe3=%le\n",&itemp,color_,&demiaxe1,&demiaxe2,&demiaxe3);    
    color = color_;
    shapeId_ = existPlanxShape(demiaxe1,demiaxe2,demiaxe3);
    
    if(shapeId_==-1)//le shape n'existe pas
    {
        shapeProp shp_pro;        
        Shape shape = createPlanxShape(demiaxe1,demiaxe2,demiaxe3);
        double density = densities_[material];
        id_shape_PLANx_.push_back(vec3r(demiaxe1,demiaxe2,demiaxe3));
        
        shp_pro.ishape = iShapes_;
        shp_pro.avrd = avrd;
        shp_pro.I1 = I1;
        shp_pro.I2 = I2;
        shp_pro.I3 = I3;
        shp_pro.density = density;
        shape_pro.push_back(shp_pro);

        Shapes.push_back(shape);
        index_shapes_.push_back(shape.name);
        iShapes_++;
    }
    else
    {
        std::string  name="PLANx" + std::to_string(shapeId_);
        index_shapes_.push_back(name);
    }
    P.homothety = 1.;
}

void Lmgc90::readPOLYR(std::ifstream & in, Particle & P, std::string & color, std::string & material, double & avrd, double & I1, double & I2, double &I3)
{
    std::string stemp,line;
    int itemp;
    int shapeId_ = -1;
    unsigned int nv=0,nf=0;
    
    Shape shape = Shape();
    shape.name = "POLYR" + std::to_string(id_shape_POLYR_);
    shape.radius = avrd/10.;
    shape.preCompDone='y';
    
    in>>itemp>>stemp>>color>>stemp>>nv>>stemp>>nf;
    std::getline(in, line);
    for(size_t iv=0;iv<nv;iv++)
    {
        vec3r vertex;
        std::getline(in, line);
        std::sscanf(line.c_str()," coo1=%le  coo2=%le  coo3=%le\n",&vertex.x,&vertex.y,&vertex.z);
        shape.vertex.push_back(vertex);
    }
    
    size_t v1,v2,v3;
    
    for (size_t e = 0; e < nf; ++e) {
        std::vector<size_t> ids;
        std::getline(in, line);
        std::sscanf(line.c_str()," ver1= %li  ver2= %li  ver3= %li\n",&v1,&v2,&v3);
        
        v1--;v2--;v3--;
                
        //il faut surement netoyer les doublons dans edge
        addEgde(shape.edge,std::pair<size_t, size_t>(v1, v2));
        addEgde(shape.edge,std::pair<size_t, size_t>(v2, v3));
        addEgde(shape.edge,std::pair<size_t, size_t>(v3, v1));
        
        ids.push_back(v1);
        ids.push_back(v2);
        ids.push_back(v3);
        shape.face.push_back(ids);
    }
    
    //adapter les vertexes a la bonne taille
    //reduceSizeShape(shape,shape.radius);

    //A utiliser uniquement lorsqu on fais es polyhedres reguliers
    int netoyage=0;
    do
    {
        netoyage=netoyagePolyr(shape);
        shape.defineVertexConnectivity();
    }
    while (netoyage!=0);
    
    shapeProp shp_pro;
    double density = densities_[material];
    shp_pro.ishape = iShapes_;
    shp_pro.avrd = avrd;
    shp_pro.I1 = I1;
    shp_pro.I2 = I2;
    shp_pro.I3 = I3;
    shp_pro.density = density;
    // ici soit j'ajoute un noveau POLYR
    // soit j affecte un shape existant a la
    // la particule
    shapeId_ = existPolyrShape(shape,shp_pro);
    
    if(shapeId_==-1)//le shape n'existe pas
    {
        shape_pro.push_back(shp_pro);
        index_shapes_.push_back(shape.name);
        avrd_polyr_.push_back(avrd);
        Shapes.push_back(shape);
        iShapes_++;
        id_shape_POLYR_++;
        P.homothety = 1.;
    }
    else
    {
        std::string  name="POLYR" + std::to_string(shapeId_);
        index_shapes_.push_back(name);
        P.homothety = avrd/avrd_polyr_[shapeId_];
    }
}

Shape Lmgc90::createSpherShape(double radius)
{
    Shape S = Shape();
    S.name="SPHER";
    S.radius=radius;
    S.preCompDone='y';
    S.volume = (4./3.)*Mth::pi*radius*radius*radius;
    S.inertia_mass = vec3r(.4*radius*radius,.4*radius*radius,.4*radius*radius);
    S.obb.extent = vec3r(radius,radius,radius);
    S.obb.center = vec3r(0.,0.,0.);
    S.obb.e[0] = vec3r(1.,0.,0.);
    S.obb.e[1] = vec3r(0.,1.,0.);
    S.obb.e[2] = vec3r(0.,0.,1.);
    vec3r pos=S.obb.center;
    S.vertex.push_back(pos);
    
    return S;
}

Shape Lmgc90::createPlanxShape(double& axe1,double& axe2,double& axe3)
{
    Shape S = Shape();
    S.name="PLANx" + std::to_string(id_shape_PLANx_.size());
    S.radius = axe3/2.;
    vec3r v[8];
    
    v[0] = vec3r(axe1,axe2,0);
    v[1] = vec3r(axe1,-axe2,0);
    v[2] = vec3r(-axe1,-axe2,0);
    v[3] = vec3r(-axe1,axe2,0);

    for (size_t i = 0; i < 4; ++i) S.vertex.push_back(v[i]);

    //edges d'en bas
    S.edge.push_back(std::pair<size_t, size_t>(0, 1)); //0
    S.edge.push_back(std::pair<size_t, size_t>(1, 2)); //1
    S.edge.push_back(std::pair<size_t, size_t>(2, 3)); //2
    S.edge.push_back(std::pair<size_t, size_t>(3, 0)); //3

    std::vector<size_t> ids0;
    ids0.push_back(0);
    ids0.push_back(1);
    ids0.push_back(2);
    ids0.push_back(3);
    S.face.push_back(ids0);

    return S;
}

/*
Shape Lmgc90::createPlanxShape(double& axe1,double& axe2,double& axe3)
{
    Shape S = Shape();
    S.name="PLANx" + std::to_string(id_shape_PLANx_.size());
    
    vec3r v[8];
    
    v[0] = vec3r(axe1,axe2,-axe3);
    v[1] = vec3r(axe1,-axe2,-axe3);
    v[2] = vec3r(-axe1,-axe2,-axe3);
    v[3] = vec3r(-axe1,axe2,-axe3);

    v[4] = vec3r(axe1,axe2,axe3);
    v[5] = vec3r(axe1,-axe2,axe3);
    v[6] = vec3r(-axe1,-axe2,axe3);
    v[7] = vec3r(-axe1,axe2,axe3);
    
    for (size_t i = 0; i < 8; ++i) S.vertex.push_back(v[i]);

    //edges d'en bas
    S.edge.push_back(std::pair<size_t, size_t>(0, 1)); //0
    S.edge.push_back(std::pair<size_t, size_t>(1, 2)); //1
    S.edge.push_back(std::pair<size_t, size_t>(2, 3)); //2
    S.edge.push_back(std::pair<size_t, size_t>(3, 0)); //3

    //edges d'en haut
    S.edge.push_back(std::pair<size_t, size_t>(4, 5)); //4
    S.edge.push_back(std::pair<size_t, size_t>(5, 6)); //5
    S.edge.push_back(std::pair<size_t, size_t>(6, 7)); //6
    S.edge.push_back(std::pair<size_t, size_t>(7, 4)); //7
    
    //edges verticales
    S.edge.push_back(std::pair<size_t, size_t>(0, 4)); //8
    S.edge.push_back(std::pair<size_t, size_t>(1, 5)); //9
    S.edge.push_back(std::pair<size_t, size_t>(2, 6)); //10
    S.edge.push_back(std::pair<size_t, size_t>(3, 7)); //11

    std::vector<size_t> ids0;
    ids0.push_back(0);
    ids0.push_back(1);
    ids0.push_back(2);
    ids0.push_back(3);
    S.face.push_back(ids0);
    S.faceArea.push_back(1.);

    std::vector<size_t> ids1;
    ids1.push_back(4);
    ids1.push_back(5);
    ids1.push_back(6);
    ids1.push_back(7);
    S.face.push_back(ids1);
    S.faceArea.push_back(1.);

    std::vector<size_t> ids2;
    ids2.push_back(0);
    ids2.push_back(1);
    ids2.push_back(5);
    ids2.push_back(4);
    S.face.push_back(ids2);
    S.faceArea.push_back(1.);

    std::vector<size_t> ids3;
    ids3.push_back(1);
    ids3.push_back(2);
    ids3.push_back(6);
    ids3.push_back(5);
    S.face.push_back(ids3);
    S.faceArea.push_back(1.);

    std::vector<size_t> ids4;
    ids4.push_back(2);
    ids4.push_back(3);
    ids4.push_back(7);
    ids4.push_back(6);
    S.face.push_back(ids4);
    S.faceArea.push_back(1.);

    std::vector<size_t> ids5;
    ids5.push_back(0);
    ids5.push_back(3);
    ids5.push_back(7);
    ids5.push_back(4);
    S.face.push_back(ids5);
    S.faceArea.push_back(1.);


    return S;
}
*/

void Lmgc90::updateshapeId()
{
    for (size_t t = 0; t < Shapes.size(); ++t) {
        shapeId[Shapes[t].name] = t;
    }
}

int Lmgc90::existPlanxShape(double &axe1,double &axe2,double &axe3)
{
    int notExist = -1;
    
    for (size_t t = 0; t < Shapes.size(); ++t)
    {
        Shape & shape = Shapes[t];
        std::string nameShape = shape.name;
        
        if(truncateLongString(nameShape,5)==PLANx)
        {
            vec3r extent = vec3r(axe1,axe2,axe3);
            int id_shape = extractnumber(shape.name);
            if(extent==id_shape_PLANx_[id_shape]) return id_shape;
        }
    }
    return notExist;
}

int Lmgc90::existPolyrShape(Shape& shape,shapeProp& shapepro)
{
    int notExist = -1;

    vec3r  I = vec3r(shapepro.I1,shapepro.I2,shapepro.I3);    
    I /= shapepro.avrd * shapepro.avrd * shapepro.avrd * 4.*M_PI/3.;;

    for (size_t t = 0; t < Shapes.size(); ++t)
    {
        Shape & shp = Shapes[t];
        shapeProp shpro = shape_pro[t];
        
        vec3r  It = vec3r(shpro.I1,shpro.I2,shpro.I3);
        It /= shpro.avrd * shpro.avrd * shpro.avrd * 4.*M_PI/3.;

        std::string nameShape = shp.name;

        if(truncateLongString(nameShape,5)==POLYR)
        {
            int id_shape = extractnumber(shp.name);

            if( (shape.vertex.size()==shp.vertex.size()) && (shape.edge.size()==shp.edge.size()) && (shape.face.size()==shp.face.size()))
            {   
                //je rajoute aussi une autre condition sur les moments d inertie
                if( (It.x/I.x) ==  (It.y/I.y) && (It.x/I.x) == (It.z/I.z) ) return id_shape;
            }
        }
    }
    return notExist;
}

unsigned int Lmgc90::extractnumber(std::string & str )
{    
    int number;
    std::string temp;
    
    for (unsigned int i=0; i < str.size(); i++)
        if (isdigit(str[i])) temp += str[i];
    
    std::istringstream stream(temp);
    stream >> number;
    return number;
}

//Lhassan
/**
 @brief Load a configuration file named name
 @param[in]  name  The name of the conf-file
 */
void Lmgc90::loadConfLMGC90(const char* name) {
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

    parser.kwMap["knContact"] = __DO__(conf){
        conf>>kn;
    };   
    parser.kwMap["en2Contact"] = __DO__(conf){
        conf>>en;
    };  
    parser.kwMap["ktContact"] = __DO__(conf){
        conf>>kt;
    };  
    parser.kwMap["krContact"] = __DO__(conf){
        conf>>kr;
    };   
    parser.kwMap["murContact"] =  __DO__(conf){
        conf>>mur;
    };   
    parser.kwMap["iconf"] = __GET__(conf, iconf);
    parser.kwMap["nDriven"] = __GET__(conf, nDriven);

    parser.kwMap["correspondingSphere"] = __DO__(conf) {
        size_t value;
        conf >> value;
        if(value==1) correspondingSphere = true;
    };

    parser.kwMap["radius"] = __GET__(conf, rInput);
    parser.kwMap["aabbextent"] = __GET__(conf, aabbextent);    

    // This single line actually parses the file
    parser.parse(conf);
    if (conf.is_open()) conf.close();
}

//LHASSAN
/**
 @brief Save a input.txt '
 */
void Lmgc90::saveConfLMGC90(){
    char fname[256];
    sprintf(fname, "input.txt");
    std::ofstream conf(fname);
    
    conf << "Rockable " << CONF_VERSION_DATE << '\n';  // format: progName version-date(dd-mm-yyyy)
    conf << "t " << t << '\n';
    conf << "tmax " << tmax << '\n';
    conf << "dt " << dt << '\n';
    conf << "interVerlet " << interVerlet << '\n';
    conf << "interConf " << interConf << '\n';
    conf << "DVerlet " << DVerlet << '\n';
    conf << "dVerlet " << dVerlet << '\n';
    for (size_t grp = 0; grp < properties.ngroup; grp++) {
        double density = properties.get(idDensity, grp);
        if (density > 0.0)  // TODO implement properties.isDefined(idDensity, grp)
            conf << "density " << grp << " " << density << '\n';
    }
    conf << "gravity " << gravity << '\n';
    
    #ifdef PERIODIC_BOX
    computeAABB();
    periodicity = aabb.max - aabb.min;
    conf << "periodicity " << periodicity << '\n';
    conf << "box_sig " << b.sig << '\n';//@LA
    conf << "box_vel " << b.vel << '\n';//@LA
    conf << "box_damp_coef " << b.coefm << '\n';//@LA
    conf << "box_mas_Totmass " << b.mass << '\n';//@LA
    #endif

    if (bodyForce != nullptr) {
        conf << "BodyForce ";
        bodyForce->write(conf);
    }

    for (auto it = optionNames.begin(); it != optionNames.end(); ++it) {
        if(it->first=="forceLaw")
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
    conf << "shapeFile " << "shape.shp" << '\n';
    if (CommBox().sep != ' ') conf << "separator " << CommBox().keywordFromSep() << '\n';
    conf << std::scientific << std::setprecision(CommBox().precision);
    
    size_t nbp = Particles.size();
    conf << "Particles " << nbp << '\n';
    
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
    
    if (conf.is_open()) conf.close();
}

void Lmgc90::reduceSizeShape(Shape & shape, double & radius)
{
    std::string str = shape.name;
    if(truncateLongString(str,5)==PLANx)
    {
        std::vector<vec3r> & v = shape.vertex;
        v[0] -= vec3r(radius,radius,0);
        v[1] -= vec3r(radius,-radius,0);
        v[2] -= vec3r(-radius,-radius,0);
        v[3] -= vec3r(-radius,radius,0);
    }
    else
    for(size_t iv=0;iv<shape.vertex.size();iv++)
    {
        vec3r & X = shape.vertex[iv];
        vec3r n = X;
        n.normalize();
        X -= n * radius;
    }
}

double volumeAABB(OBB & aa)
{
    return 2.*(aa.extent.x)*(aa.extent.y)*(aa.extent.z); 
}

double volumeAABB(AABB & aa)
{
    vec3r vv = (aa.max - aa.min);
    return (vv.x)*(vv.y)*(vv.z); 
}

void Lmgc90::exacteSize(Particle & P)
{
    // Ici, j'essaye de trouver le volume exacte
    double epsilon = 1.e-13; 

    //on cree un shape aux qui va etre modifie
    Particle aux = P;

    double maxH = 5.;
    double minH = 0.5;
    double H = 0.;

    //La boite qui contient le shape initial avec V0 le volume initial
    P.updateObb();
    OBB box = P.obb;
    double V0 = volumeAABB(box);

    Hmax_[&P] = P.homothety;

    aux.homothety = minH;
    aux.updateObb();
    OBB ii = aux.obb;
    ii.enlarge(aux.MinskowskiRadius());

    double V = volumeAABB(ii);
    
    long maxIteration =50000;
    size_t i =0;

    while (fabs(V0-V) > epsilon && (long)i<maxIteration )
    {
        H = (maxH+minH)/2.;
        aux.homothety = H;
        aux.updateObb();
        OBB mm = aux.obb;
        mm.enlarge(aux.MinskowskiRadius());
        V = volumeAABB(mm);

        if(V>V0)
        {
            maxH = H;
        }
        else
        {
            minH = H;
        }
        
        i++;
        aux.homothety = 1./H;
    }

    std::cout<<"shap name "<<P.shape->name<<" |V0-V| "<<fabs(volumeAABB(box)-V)<<", epsilon "<<epsilon<<std::endl;

    P.homothety = H;
}


void Lmgc90::exacteSize(Shape & shape)
{
    std::string str = shape.name;
    double radius = shape.radius;
    if(truncateLongString(str,5)==PLANx)
    {
        std::vector<vec3r> & v = shape.vertex;
        v[0] -= vec3r(radius,radius,0);
        v[1] -= vec3r(radius,-radius,0);
        v[2] -= vec3r(-radius,-radius,0);
        v[3] -= vec3r(-radius,radius,0);
    }
    else
    {
        // Ici, j'essaye de trouver le volume exacte
        double epsilon = 1.e-16; 
        //on cree un shape aux qui va etre modifie
        
        double maxH = 5.;
        double minH = 0.5;
        double H = 0.;
        //La boite qui contient le shape initial avec V0 le volume initial

        AABB box(shape.vertex); 

        std::cout<<"------------------------------------------------------------------"<<std::endl<<std::endl;
        std::cout<<"Shape name : "<<shape.name<<std::endl;
        vec3r originextent = (box.max-box.min)/2.;
        std::cout<<"original obb extent : "<<originextent.x<<" "<<originextent.y<<" "<<originextent.z<<std::endl;


        double V0 = volumeAABB(box);

        shape.homothety(minH);
        AABB ii(shape.vertex);
        ii.enlarge(shape.radius);
        double V = volumeAABB(ii);
    
        long maxIteration =50000;
        size_t i =0;

        while (fabs(V0-V)/V0 > epsilon && (long)i<maxIteration )
        {
            H = (maxH+minH)/2.;
            shape.homothety(H);
            
            AABB mm(shape.vertex);
            mm.enlarge(shape.radius);
            V = volumeAABB(mm);

            if(V>V0)
            {
                maxH = H;
            }
            else
            {
                minH = H;
            }
        
            i++;
            shape.homothety(1./H);
        }

        std::cout<<" |V0-V| "<<fabs(V0-V)/V0<<" / "<<epsilon<<" homothety = "<<H << std::endl;
        shape.homothety(H);

        AABB finalbox(shape.vertex);
        finalbox.enlarge(radius);
        shape.obb.extent = (finalbox.max - finalbox.min)/2.;
        std::cout<<"final obb extent : "<<shape.obb.extent.x<<" "<<shape.obb.extent.y<<" "<<shape.obb.extent.z<<std::endl;

        std::cout<<"------------------------------------------------------------------"<<std::endl;

    }
}

int Lmgc90::netoyagePolyr(Shape &shape)
{
    double epsilon = 1.e-13; //epsilon pour le produit scalaire des deux normales
    double dotscalar = 0.;   //au plan
    
    for(size_t i=0;i<shape.face.size();i++)
    {
        vec3r ni  =  get_normal_face(i,shape);
         
        for(size_t j=i+1;j<shape.face.size();j++)
        {
            if(siTowFacesConnected(shape.face[i],shape.face[j]))
            {
                vec3r nj = get_normal_face(j,shape);                
                dotscalar = fabs(fabs(dot(ni,nj))-1.);
                if(dotscalar<=epsilon)
                {
                    //là il faut penser à une methode de fusion des faces et enlever l'edge commun
                    
                    //decider quel face il faut suprimer d abord
                    //ici c'est la face qui est en deuxieme parametre
                    std::pair<size_t,size_t> edge;
                    
                    if(shape.face[i].size()>shape.face[j].size())
                        edge = commonEdgeFacesConnected(shape.face[i],shape.face[j]);
                    else
                        edge = commonEdgeFacesConnected(shape.face[j],shape.face[i]);
                    
                    removeEgde(shape.edge,edge);
                    
                    //concatener les deux faces
                    
                    if(shape.face[i].size()>shape.face[j].size())
                    {
                        concatener_deux_faces(shape.face,i,j,edge);   
                    }
                    else
                    {
                        concatener_deux_faces(shape.face,j,i,edge);
                    }
                    return 1;
                }
            }
        }
    }
    return 0;
}

bool Lmgc90::siTowFacesConnected(std::vector<size_t> face1,std::vector<size_t> face2)
{
    size_t nb_common_vertex = 0;
    
    for(size_t i=0;i<face1.size();i++)
    {
        for(size_t j=0;j<face2.size();j++)
        {
            if(face1[i]==face2[j]) nb_common_vertex++;
        }
    }
    
    if(nb_common_vertex>=2) return true;
    return false;
}

std::pair<size_t,size_t> Lmgc90::commonEdgeFacesConnected(std::vector<size_t> face1,std::vector<size_t> face2)
{
    std::pair<size_t,size_t> common_edge;
    size_t nb_common_vertex = 0;
    
    for(size_t i=0;i<face1.size();i++)
    {
        for(size_t j=0;j<face2.size();j++)
        {
            if(face1[i]==face2[j])
            {
                nb_common_vertex++;
                if(nb_common_vertex==1)common_edge.first = face2[j];
                if(nb_common_vertex==2)common_edge.second = face2[j];
            }
        }
    }
    return common_edge;
}

void Lmgc90::removeEgde(std::vector<std::pair<size_t, size_t> > &vedge, std::pair<size_t, size_t> edge)
{
    std::vector<std::pair<size_t, size_t> >::iterator it = vedge.begin();
    while(it != vedge.end())
    {
        std::pair<size_t, size_t> edg = (std::pair<size_t, size_t>)*it;
        if( ((edg.first==edge.first) && (edg.second==edge.second)) || ((edg.first==edge.second) && (edg.second==edge.first))  )
        {
            vedge.erase(it);
            break;
        }
        else ++it;
    }
}

void Lmgc90::removeFace(std::vector<std::vector<size_t> > & face,size_t f1)
{
    // erase the f1 th element
    face.erase (face.begin()+f1);
}

void Lmgc90::addEgde(std::vector<std::pair<size_t, size_t> > &vedge, std::pair<size_t, size_t> edge)
{
    for(size_t i =0;i<vedge.size();i++)
    {
        std::pair<size_t, size_t> edg = vedge[i];
        if( ((edg.first==edge.first) && (edg.second==edge.second)) || ((edg.first==edge.second) && (edg.second==edge.first))  ) return;
    }
    vedge.push_back(edge);
}

void Lmgc90::concatener_deux_faces(std::vector<std::vector<size_t> > & face,size_t fmax,size_t fmin, std::pair<size_t,size_t> edge)
{
    //edge contient l'orientation de face max
    
    int debut=0;
    int  fin=0;
    //il me faut trouver la position de edge.first dans fmin
    for(size_t i =0;i<face[fmin].size();i++)
    {
        if(edge.first==face[fmin][i]) debut = i;
        if(edge.second==face[fmin][i]) fin  = i;
    }
    
    std::vector<int> append;
    
    int nbvmax = face[fmin].size()-2;//nombre de vertex maximal a ajouter
    
    append.resize(nbvmax);
    
        if(fin-debut>0) //le cas edial les deux se suivent
        {
            if(fin-debut==1)
            {
                int nbv = 0;
                for(int i=debut-1;i>=0;i--)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            
                for(int i=face[fmin].size()-1;i>fin;i--)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            }
            else //extrimite
            {
                int nbv = 0;
                for(int i=debut+1;i<fin;i++)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            }
        }
        else //le cas les deux se suivnet mais il faut inverser
        {
            if(fin-debut==-1)
            {
                int nbv = 0;
                for(size_t i=debut+1;i<face[fmin].size();i++)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            
                for(int i=0;i<fin;i++)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            }
            else //extrimite
            {
                int nbv = 0;
                for(int i=debut-1;i>fin;i--)
                {
                    append[nbv] = face[fmin][i];
                    nbv++;
                }
            }
        }
    
    std::vector<int> ids;
    
    //il me faut trouver la position de edge.first dans fmin
    for(size_t i =0;i<face[fmax].size();i++)
    {
        if(edge.first==face[fmax][i]) debut = i;
        if(edge.second==face[fmax][i]) fin  = i;
    }
    
    if(fin-debut==1)
    {
        //les deux se suivent
        for(int i =0;i<=debut;i++)
        {
            ids.push_back(face[fmax][i]);
        }
        for(size_t i =0;i<append.size();i++)
        {
            ids.push_back(append[i]);
        }
        for(size_t i =fin;i<face[fmax].size();i++)
        {
            ids.push_back(face[fmax][i]);
        }
    }
    else
    {
        //les deux sont à l'extremmité
        for(size_t i =0;i<face[fmax].size();i++)
        {
            ids.push_back(face[fmax][i]);
        }
        for(size_t i =0;i<append.size();i++)
        {
            ids.push_back(append[i]);
        }
    }
    
    face[fmax].clear();
    for(size_t i =0;i<ids.size();i++)
    {
        face[fmax].push_back(ids[i]);
    }
    
    //j enleve la face
    removeFace(face,fmin);
}

vec3r Lmgc90::get_normal_face( size_t & f, Shape & shape)
{
    vec3r n;
    std::vector<size_t> face = shape.face[f];

    //je cherche les deux vecteur qui ont le plus grand produit scalaire :
    std::vector<vec3r> vecteurs; 

    for(size_t i=0;i<face.size();i++)
    {
        vec3r vi = shape.vertex[face[i]];
        for(size_t j=i+1;j<face.size();j++)
        {
            vec3r vj = shape.vertex[face[j]];
            vec3r vji = vj - vi;
            vecteurs.push_back(vji);
        }
    }

    double maxdot = 0.;
    size_t maxi =0,maxj=0;

    for(size_t i=0;i<vecteurs.size();i++)
    {
        vec3r vij_i = vecteurs[i];
        for(size_t j=i+1;j<vecteurs.size();j++)
        {
            vec3r vij_j = vecteurs[j];

            if(fabs(dot(vij_i,vij_j))>maxdot)
            {
                maxdot = fabs(dot(vij_i,vij_j));
                maxi = i;
                maxj = j;
            }
        }
    }

    n = vecteurs[maxi] ^ vecteurs[maxj];
    n.normalized();

    return n; 
}

void Lmgc90::ToDEM3DLiens()
{
    //je vais ecrire chaque lien origine ext

    std::string fout_ = "./Spl.dat";
    std::ofstream out(fout_.c_str());
    
    if(out)
    {
        out << Particles.size()<<std::endl;

        for(size_t p = 0; p < Particles.size(); ++p)
        {
            vec3r orig = Particles[p].pos;
            double Rj = Particles[p].MinskowskiRadius();
            
            out << Particles[p].shape->vertex.size() <<" "<<Rj<<" "<< orig.x<<" "<<orig.y<<" "<<orig.z<<std::endl;     
            for (size_t iv = 0; iv < Particles[p].shape->vertex.size(); iv++) 
            {
                double h1 = Particles[p].homothety;
                vec3r ext1 = Particles[p].GlobVertex(iv);
                Particles[p].homothety = Hmax_[&Particles[p]];    
                vec3r ext2 = Particles[p].GlobVertex(iv);
                double hmax = norm(ext2-orig);    
                Particles[p].homothety = h1;
                out <<hmax<<" "<<ext1.x<<" "<<ext1.y<<" "<<ext1.z<<std::endl; 
            }
        } 
        out.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fout_.c_str() <<  std::endl;
    }
}

/*
//Ajustement pour les PArticules
void Lmgc90::ToDEM3DFaces()
{
    //je vais ecrire chaque lien origine ext

    std::string fout_ = "./Spl.dat";
    std::ofstream out(fout_.c_str());
    
    
    msg::bestPrecision(out);

    if(out)
    {
        out << Particles.size()<<std::endl;

        for(size_t p = 0; p < Particles.size(); ++p)
        {
            vec3r orig = Particles[p].pos;
            double Rj = Particles[p].MinskowskiRadius();
            
            out << Particles[p].shape->vertex.size() <<" "<<Rj<<" "<< orig.x<<" "<<orig.y<<" "<<orig.z<<std::endl;     
            for (size_t iv = 0; iv < Particles[p].shape->vertex.size(); iv++) 
            {
                double h1 = Particles[p].homothety;
                vec3r ext1 = Particles[p].GlobVertex(iv);
                Particles[p].homothety = Hmax_[&Particles[p]];    
                vec3r ext2 = Particles[p].GlobVertex(iv);
                double hmax = norm(ext2-orig);    
                Particles[p].homothety = h1;
                out <<hmax<<" "<<ext1.x<<" "<<ext1.y<<" "<<ext1.z<<std::endl; 
            }

            out << Particles[p].shape->face.size()<<std::endl;     
            for (size_t f = 0; f < Particles[p].shape->face.size(); f++) 
            {
                std::vector<size_t> face = Particles[p].shape->face[f];
                out <<face.size();
                for(size_t fi = 0; fi < face.size(); fi++){out<<" "<<face[fi];}
                out<<std::endl;
            }
        } 
        out.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fout_.c_str() <<  std::endl;
    }
}
*/

//Ajustement pour les Shapes (juste pour les polyres)
void Lmgc90::ToDEM3DFaces()
{
    //je vais ecrire chaque lien origine ext 
    std::string fout_ = "./Spl_in.dat";
    std::ofstream out(fout_.c_str());
    
    msg::bestPrecision(out);

    if(out)
    {
        out << id_shape_POLYR_ <<std::endl;

        for(size_t p = 0; p < Shapes.size(); ++p)
        {
            std::string str = Shapes[p].name;
            if((truncateLongString(str,5)==POLYR))    
            {
                vec3r orig = Shapes[p].position;
                double Rj = Shapes[p].radius;
                out << Shapes[p].name << std::endl; 
                out << Shapes[p].vertex.size() <<" "<<Rj<<" "<< orig.x<<" "<<orig.y<<" "<<orig.z<<std::endl;     
                for (size_t iv = 0; iv < Shapes[p].vertex.size(); iv++) 
                {
                    vec3r ext = Shapes[p].vertex[iv];
                    double hmax = norm(ext-orig);    
                    out <<hmax<<" "<<ext.x<<" "<<ext.y<<" "<<ext.z<<std::endl; 
                }

                out << Shapes[p].face.size()<<std::endl;     
                for (size_t f = 0; f < Shapes[p].face.size(); f++) 
                {
                    std::vector<size_t> face = Shapes[p].face[f];
                    out <<face.size();
                    for(size_t fi = 0; fi < face.size(); fi++){out<<" "<<face[fi];}
                    out<<std::endl;
                }
            }
        } 
        out.close();
    }
    else
    {
        std::cerr << "Error message : " << "Impossible de lire le fichier : " << fout_.c_str() <<  std::endl;
    }
}

size_t Lmgc90::shape_by_name(std::string & name)
{
    for(size_t p = 0; p < Shapes.size(); ++p)
    {
        if(Shapes[p].name == name) return p;
    }
    return -1;
}


//Ajustement pour les Shapes
void Lmgc90::FromDEM3DFaces()
{
    //je vais ecrire chaque lien origine ext
    std::string line;
    std::string fin_ = "./Spl_out.dat";
    std::ifstream fin(fin_.c_str());
    
    double r;
    vec3r  o;
    size_t n;
    
    size_t np;

    if(fin)
    {
        fin >> np;
        
        for(size_t p = 0; p < np; ++p)
        {
            std::string str;
            fin>>str;
            size_t i = shape_by_name(str);
            Shape & shp = Shapes[i];
            fin >> n >> r>> o.x >>  o.y >> o.z; 
            std::cout<<n<<" "<<r<<" "<<o.x<<" "<<o.y<<" "<<o.z<<std::endl;  
            for (size_t iv = 0; iv < shp.vertex.size(); iv++) 
            {
                vec3r & ext = shp.vertex[iv];
                fin >> ext.x>>ext.y>>ext.z; 
                shp.radius = r;
            }
        } 
        fin.close();
    }
}