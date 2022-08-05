/**
 * @file
 * @brief Implementation of the VtkPolytopicReader
 * @author Tarzis Maurer
 * @date   2022-04-12
 * @copyright ETH Zurich
 */

#include "vtk_polytopic_reader.h"

#include <iostream>
#include <fstream>


namespace lf::io {

std::vector<std::string> splitString(std::string s){
    std::vector<std::string> v;
	
	std::string temp = "";
	for(int i=0;i<s.length();++i){
		
		if(s[i]==' '){
			v.push_back(temp);
			temp = "";
		}
		else{
			temp.push_back(s[i]);
		}
		
	}
	v.push_back(temp);
	return v;
}

VtkPolytopicReader::VtkPolytopicReader(std::unique_ptr<mesh::polytopic2d::MeshFactory> factory,
                    const std::string& filename) : mesh_factory_(std::move(factory)) {
    
    using coord_t = Eigen::Vector2d;

    //Note: currently only implemented in 2D
    dim_t dim_mesh = mesh_factory_->DimMesh();
    dim_t dim_world = mesh_factory_->DimWorld();

    LF_VERIFY_MSG(
      dim_mesh == 2 && dim_world == 2,
      "VtkPolytopicReader currentl supports only 2D.");

    //preparation to read from file
    std::ifstream inFile;
    inFile.open(filename, std::ios::in);
    LF_VERIFY_MSG(inFile, "Unable to open file");

    //section which is currently read
    VtkSection section = header;

    //placeholders fopr strings that are read:
    //the whole line and the line split into words separated by whitespace
    std::string line;
    std::vector<std::string> split;

    //number of points in the mesh
    size_type number_points;
    //number of polygons in the mesh
    size_type number_polygons;
    //number after number polygons
    size_type number_point_index;
    //vector of offset numbers
    std::vector<size_type> offset_vector;
    //vector of connectivity info
    std::vector<size_type> connectivity_vector;

    //read the file
    while(true){

        switch(section){
            case(header):
                std::getline(inFile, line);
                split = splitString(line);
                if (split[0] == "POINTS"){
                    number_points = std::stoul(split[1]);
                    section = points;
                    break;
                }
                break;

            case(points):
                std::getline(inFile, line);
                if (line.size() == 1){ // empty line => section is finished
                    section = metadata;
                    break;
                }
                split = splitString(line);
                if (split[0] == "POLYGONS"){
                    number_polygons = std::stoul(split[1]);
                    number_point_index = std::stoul(split[2]);
                    section = offsets;
                    offset_vector.reserve(number_polygons);
                    //std::getline(inFile, line); // next line is not used
                    break;
                }
                //read first two coordinates of each point and add the point to the MeshFactory
                for (int i = 0; i < split.size() - 1 ; i += 3){
                    double x = std::stod(split[i]);
                    double y = std::stod(split[i+1]);
                    mesh_factory_->AddPoint(coord_t({x, y}));
                    LF_VERIFY_MSG(std::stod(split[i+2]) == 0.0, "Something wrong with order of coordinates, z coordinate should be 0.0");
                }
                break;

            case(metadata):
                std::getline(inFile, line);
                split = splitString(line);
                if (split[0] == "POLYGONS"){
                    number_polygons = std::stoul(split[1]);
                    number_point_index = std::stoul(split[2]);
                    section = offsets;
                    offset_vector.reserve(number_polygons);
                    std::getline(inFile, line); // next line is not used
                    break;
                }
                break;

            case(offsets):
                std::getline(inFile, line); 
                split = splitString(line);
                if (split[0] == "CONNECTIVITY"){ // section is finished
                    section = connectivity;
                    break;
                } else if (split[0] == "OFFSETS"){ //line is not used
                    break;
                }
                for (int i = 0; i < split.size() - 1; i++){       //last string is no number
                    offset_vector.push_back(std::stoul(split[i]));
                }
                break;

            case(connectivity):
                std::getline(inFile, line); 
                split = splitString(line);
                if (line.size() == 1){ // empty line => section is finished
                    goto exit_switch;
                }
                for (int i = 0; i < split.size() - 1; i++){       //last string is no number
                    connectivity_vector.push_back(std::stoul(split[i]));
                }
                break;
                
        }

        exit_switch: ;
        if(inFile.eof()){
            break;
        }
    }

    inFile.close();

    //now construct the polygons
    //loop over all polygons
    for (int i = 0; i < offset_vector.size() - 1; i++){
        int start = offset_vector.at(i);
        int end = offset_vector.at(i+1);
        const int number_nodes = end - start;
        std::vector<size_type> node_indices(number_nodes);
        int counter = 0;
        //loop over indices of the polygon
        for(start; start < end; start++){
            size_type next_node = connectivity_vector.at(start);
            node_indices.at(counter) = next_node;
            counter++; 
        }
        //add the polygon
        mesh_factory_->AddEntity(lf::base::RefEl::kPolygon(), node_indices, nullptr);
    }
    mesh_ = mesh_factory_->Build();
}


} // nemaspace lf::io