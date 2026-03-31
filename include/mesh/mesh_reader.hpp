#pragma once

#include <string>

class MeshReader
{
public:
    MeshReader(const std::string& filename);
    void read();


private : 
    std::string filename_;
};