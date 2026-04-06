#include "mesh/mesh_compute.hpp"
#include <cmath>
#include <unordered_map>


MeshCompute::MeshCompute(const MeshReader& reader)
    : reader_(reader)
{
    std::cout << "MeshCompute created " << std::endl;
}

void MeshCompute::compute_center_cells()
{
    const auto& X = reader_.get_CoordX();
    const auto& Y = reader_.get_CoordY();
    const auto& Z = reader_.get_CoordZ();
    const auto& T = reader_.get_Triangles();

    centercells_X.resize(T.size());
    centercells_Y.resize(T.size());
    centercells_Z.resize(T.size());

    for(size_t i=0; i< T.size(); i++)
    {
        int n0 = T[i][0];
        int n1 = T[i][1];
        int n2 = T[i][2];
        centercells_X[i] =(X[n0] + X[n1] + X[n2])/3.0;
        centercells_Y[i] =(Y[n0] + Y[n1] + Y[n2])/3.0;
        centercells_Z[i] =(Z[n0] + Z[n1] + Z[n2])/3.0;
    }

    std::cout << "Center Cells Calculated " << std::endl;
}

void MeshCompute::compute_cells_area()
{
    const auto& X = reader_.get_CoordX();
    const auto& Y = reader_.get_CoordY();
    const auto& Z = reader_.get_CoordZ();
    const auto& T = reader_.get_Triangles();

    cellsarea.resize(T.size());

    for(size_t i=0; i<T.size(); i++)
    {
        int n0 = T[i][0];
        int n1 = T[i][1];
        int n2 = T[i][2];

        double ux = X[n1] - X[n0];
        double uy = Y[n1] - Y[n0];
        double uz = Z[n1] - Z[n0];

        double vx = X[n2] - X[n0];
        double vy = Y[n2] - Y[n0];
        double vz = Z[n2] - Z[n0];

        double cx = uy * vz - uz * vy;
        double cy = uz * vx - ux * vz;
        double cz = ux * vy - uy * vx;

        cellsarea[i] = 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
    }

    std::cout << "Cells area Calculated " << std::endl;
}

void MeshCompute::compute_faces()
{
    const auto& X = reader_.get_CoordX();
    const auto& Y = reader_.get_CoordY();
    const auto& T = reader_.get_Triangles();

    faces_nodes.clear();
    faces_cells.clear();
    faces_center_x.clear();
    faces_center_y.clear();
    faces_length.clear();
    faces_normal_x.clear();
    faces_normal_y.clear();

    std::unordered_map<long long, int> face_map;

    for (int cell = 0; cell < static_cast<int>(T.size()); cell++) 
    {
        int n0 = T[cell][0];
        int n1 = T[cell][1];
        int n2 = T[cell][2];

        // Each loop we build all faces possible for the triangle
        std::array<std::array<int,2>,3> cell_faces = {{
            {n0, n1},
            {n1, n2},
            {n2, n0}
        }};

        // We sort each face to avoid dupes
        for (int k = 0; k < 3; k++)
        {
            int a = cell_faces[k][0];
            int b = cell_faces[k][1];

            int i, j;

            if (a < b)
            {
                i = a;
                j = b;
            }
            else
            {
                i = b;
                j = a;
            }

            long long key = ((long long)i << 32) | (unsigned int)j;

            auto it = face_map.find(key);

            if (it == face_map.end())
            {
                // if the face is not already in the map we add it 
                int face_id = static_cast<int>(faces_nodes.size());

                face_map[key] = face_id;

                faces_nodes.push_back({a, b});
                faces_cells.push_back({cell, -1});

                double xm = 0.5 * (X[a] + X[b]);
                double ym = 0.5 * (Y[a] + Y[b]);

                // We add the face center 
                faces_center_x.push_back(xm);
                faces_center_y.push_back(ym);

                double dx = X[b] - X[a];
                double dy = Y[b] - Y[a];
                // Same for the length
                double len = std::sqrt(dx*dx + dy*dy);
                faces_length.push_back(len);

                // Raw face normal associated with edge vector (dx, dy): (dy, -dx)
                if (len > 1e-14)
                {
                    faces_normal_x.push_back(dy/len);
                    faces_normal_y.push_back(-dx/len);
                }
                else 
                {
                    faces_normal_x.push_back(0.0);
                    faces_normal_y.push_back(0.0);
                   
                }
                
            }
            else
            {
                // if the face is already in the map we add a right cell only
                int face_id = it->second;
                faces_cells[face_id][1] = cell;
            }
        }
    }

    std::cout << "Faces computed: " << faces_nodes.size() << std::endl;
}

void MeshCompute::orient_faces()
{
    for (size_t f = 0; f < faces_nodes.size(); f++)
    {
        int cl = faces_cells[f][0];

        double nx = faces_normal_x[f];
        double ny = faces_normal_y[f];

        double vx = faces_center_x[f] - centercells_X[cl];
        double vy = faces_center_y[f] - centercells_Y[cl];

        double dot = nx * vx + ny * vy;

        if (dot < 0.0)
        {
            faces_normal_x[f] = -nx;
            faces_normal_y[f] = -ny;
        }
    }

    int bad = 0;

    for (size_t f = 0; f < faces_nodes.size(); f++)
    {
        int cl = faces_cells[f][0];

        double vx = faces_center_x[f] - centercells_X[cl];
        double vy = faces_center_y[f] - centercells_Y[cl];

        // Dot product between the face normal and the vector from the left cell center to the face center
        double dot = vx * faces_normal_x[f] + vy * faces_normal_y[f];

        if (dot <= 0.0)
        {
            bad++;
        }
        
    }

    std::cout << "Bad oriented faces = " << bad << std::endl;

    std::cout << "Faces oriented" << std::endl;
    
}

