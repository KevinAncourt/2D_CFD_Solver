#pragma once

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>


class MeshReader
{
public:
    MeshReader(const std::string& filename);
    void read_node(std::ifstream& file);
    void read_element(std::ifstream& file);
    void read_physical_names(std::ifstream& file);
    void read_entities(std::ifstream& file);
    void write_vtk(const std::string& filename) const;
    void write_vtk_wall(const std::string& filename) const;
    void write_vtk_farfield(const std::string& filename) const;
    void read();
    
    const std::vector<double>& get_CoordX() const;
    const std::vector<double>& get_CoordY() const;
    const std::vector<double>& get_CoordZ() const;
    const std::vector<std::array<int, 3>>& get_Triangles() const;
    const std::vector<std::array<int, 2>>& get_BdryEdges_farfield() const;
    const std::vector<std::array<int, 2>>& get_BdryEdges_wall() const;


private : 
    std::string filename_;
    std::vector<double> CoordX;
    std::vector<double> CoordY;
    std::vector<double> CoordZ;
    std::vector<int> lntogn;
    std::vector<int> gntoln;
    std::vector<std::array<int, 2>> edgebcwall;
    std::vector<std::array<int, 2>> edgebcfarfield;
    std::vector<std::array<int, 3>> triangles;

    void write_vtk_edges(const std::string& filename,const std::vector<std::array<int,2>>& edges) const;

    std::unordered_map<int, std::string> physicalNames1D;          // physicalTag -> "Wall"/"Farfield"
    std::unordered_map<int, std::vector<int>> curveToPhysicalTags; // entityTag(courbe) -> liste de physicalTag

    std::string get_curve_bc_name(int curveEntityTag) const;
};