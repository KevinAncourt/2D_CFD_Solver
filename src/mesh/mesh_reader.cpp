#include "mesh/mesh_reader.hpp"
#include <iostream>
#include <fstream>


MeshReader::MeshReader(const std::string& filename)
{
    filename_ = filename;
}

void MeshReader::read()
{
    std::ifstream file(filename_);  

    if (!file)
    {
        std::cerr << "Erreur : impossible d'ouvrir " << filename_ << std::endl;
        return;
    }

    std::cout << "Fichier ouvert : " << filename_ << std::endl;

    std::string line;

    // while (std::getline(file, line))
    // {
    //     std::cout << line << std::endl;  // pour tester l'affichage des lignes
    // }
}