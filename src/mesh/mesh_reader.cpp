#include "mesh/mesh_reader.hpp"


MeshReader::MeshReader(const std::string& filename)
{
    filename_ = filename;
}

void MeshReader::read_node(std::ifstream& file)
{
    std::string line;
    
    std::getline(file, line);

    std::istringstream iss(line);

    int numEntityBlocks;
    int numNodes;
    int minNodeTag;
    int maxNodeTag;

    iss >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;

    std::cout << "numEntityBlocks = " << numEntityBlocks << std::endl;
    std::cout << "numNodes        = " << numNodes << std::endl;
    std::cout << "minNodeTag      = " << minNodeTag << std::endl;
    std::cout << "maxNodeTag      = " << maxNodeTag << std::endl;

    CoordX.resize(numNodes);
    CoordY.resize(numNodes);
    CoordZ.resize(numNodes);
    lntogn.resize(numNodes);
    gntoln.resize(maxNodeTag + 1, -1);

    int localIndex = 0;

    for (int b = 0; b < numEntityBlocks; b++)
    {
        std::getline(file, line);
        std::istringstream blockHeader(line); // Convert the line into an input stream to extract values easily
 
        int entityDim;
        int entityTag;
        int parametric;
        int numNodesInBlock;

        blockHeader >> entityDim >> entityTag >> parametric >> numNodesInBlock;

        std::vector<int> nodeTags(numNodesInBlock);

        for (int i = 0; i < numNodesInBlock; i++)
        {
            std::getline(file, line);
            nodeTags[i] = std::stoi(line); // Convert the string line into an integer
        }

        for (int i = 0; i < numNodesInBlock; i++)
        {
            std::getline(file, line);
            std::istringstream coordStream(line);

            double x, y, z;
            coordStream >> x >> y >> z;

            int nodeTag = nodeTags[i];

            CoordX[localIndex] = x;
            CoordY[localIndex] = y;
            CoordZ[localIndex] = z;

            lntogn[localIndex] = nodeTag;
            gntoln[nodeTag] = localIndex;

            localIndex++;
        }
    }

    std::getline(file, line); // read $EndNodes

    std::cout << "Nodes finished" << std::endl;
                                
}

void MeshReader::read_element(std::ifstream& file)
{
    std::string line;

    std::getline(file, line);
    std::istringstream iss(line);

    int numEntityBlocks;
    int numElements;
    int minElementTag;
    int maxElementTag;

    iss >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;

    for (int iblock = 0; iblock < numEntityBlocks; iblock++)
    {
        std::getline(file, line);
        std::istringstream blockHeader(line);

        int entityDim;
        int entityTag;
        int elementType;
        int numElementsInBlock;

        blockHeader >> entityDim >> entityTag >> elementType >> numElementsInBlock;

        if (entityDim == 1 && elementType == 1)
        {
            const std::string bcName = get_curve_bc_name(entityTag);

            for (int ielem = 0; ielem < numElementsInBlock; ielem++)
            {
                std::getline(file, line);
                std::istringstream elemStream(line);

                int elementTag;
                int n1_gmsh, n2_gmsh;
                elemStream >> elementTag >> n1_gmsh >> n2_gmsh;

                int n1 = gntoln[n1_gmsh];
                int n2 = gntoln[n2_gmsh];

                if (bcName == "Wall" || bcName == "wall")
                {
                    edgebcwall.push_back({n1, n2});
                }
                else if (bcName == "Farfield" || bcName == "farfield")
                {
                    edgebcfarfield.push_back({n1, n2});
                }
                else
                {
                    std::cerr << "Warning: unknown 1D physical group for curve entityTag = "
                              << entityTag << " (name = \"" << bcName << "\")\n";
                }
            }
        }
        else if (entityDim == 2 && elementType == 2)
        {
            for (int ielem = 0; ielem < numElementsInBlock; ielem++)
            {
                std::getline(file, line);
                std::istringstream elemStream(line);

                int elementTag;
                int n1_gmsh, n2_gmsh, n3_gmsh;
                elemStream >> elementTag >> n1_gmsh >> n2_gmsh >> n3_gmsh;

                int n1 = gntoln[n1_gmsh];
                int n2 = gntoln[n2_gmsh];
                int n3 = gntoln[n3_gmsh];

                triangles.push_back({n1, n2, n3});
            }
        }
        else
        {
            for (int ielem = 0; ielem < numElementsInBlock; ielem++)
                std::getline(file, line);
        }
    }

    std::getline(file, line); // $EndElements
    std::cout << "Elements finished" << std::endl;
}

void MeshReader::write_vtk(const std::string& filename) const
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    int numPoints = CoordX.size();
    int numTriangles = triangles.size();

    out << "# vtk DataFile Version 2.0\n";
    out << "Mesh exported from MeshReader\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Points
    out << "POINTS " << numPoints << " float\n";
    for (int i = 0; i < numPoints; i++)
    {
        out << CoordX[i] << " " << CoordY[i] << " " << CoordZ[i] << "\n";
    }

    // Cells
    out << "\nCELLS " << numTriangles << " " << 4 * numTriangles << "\n";
    for (int i = 0; i < numTriangles; i++)
    {
        out << "3 "
            << triangles[i][0] << " "
            << triangles[i][1] << " "
            << triangles[i][2] << "\n";
    }

    // Cell types
    out << "\nCELL_TYPES " << numTriangles << "\n";
    for (int i = 0; i < numTriangles; i++)
    {
        out << "5\n"; // VTK_TRIANGLE
    }

    out.close();

    std::cout << "VTK file written: " << filename << std::endl;
}

void MeshReader::write_vtk_edges(const std::string& filename,
                                 const std::vector<std::array<int,2>>& edges) const
{
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Error: cannot open " << filename << std::endl;
        return;
    }

    int numPoints = CoordX.size();
    int numEdges  = edges.size();

    out << "# vtk DataFile Version 2.0\n";
    out << "Edges\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // Points
    out << "POINTS " << numPoints << " float\n";
    for (int i = 0; i < numPoints; i++)
    {
        out << CoordX[i] << " " << CoordY[i] << " " << CoordZ[i] << "\n";
    }

    // Cells
    out << "\nCELLS " << numEdges << " " << 3 * numEdges << "\n";
    for (int i = 0; i < numEdges; i++)
    {
        out << "2 "
            << edges[i][0] << " "
            << edges[i][1] << "\n";
    }

    // Types
    out << "\nCELL_TYPES " << numEdges << "\n";
    for (int i = 0; i < numEdges; i++)
    {
        out << "3\n"; // VTK_LINE
    }

    out.close();

    std::cout << "VTK edges written: " << filename << std::endl;
}

void MeshReader::write_vtk_wall(const std::string& filename) const
{
    write_vtk_edges(filename, edgebcwall);
}

void MeshReader::write_vtk_farfield(const std::string& filename) const
{
    write_vtk_edges(filename, edgebcfarfield);
}

void MeshReader::read()
{
    std::ifstream file(filename_);  

    if (!file)
    {
        std::cerr << "Error : Can't open Mesh " << filename_ << std::endl;
        return;
    }
    std::string line;

    std::cout << "Mesh Open : " << filename_ << std::endl;
    while (std::getline(file, line))
    {
        if (line == "$PhysicalNames")
        {
            read_physical_names(file);
        }
        
        else if (line == "$Entities")
        {
            read_entities(file);
        }

        if (line == "$Nodes")
        {
        MeshReader::read_node(file);
        }

        if (line == "$Elements")
        {
        MeshReader::read_element(file);
        }
        
    }
}

const std::vector<double>& MeshReader::get_CoordX() const
{
    return CoordX;
}

const std::vector<double>& MeshReader::get_CoordY() const
{
    return CoordY;
}

const std::vector<double>& MeshReader::get_CoordZ() const
{
    return CoordZ;
}

const std::vector<std::array<int, 3>>& MeshReader::get_Triangles() const
{
    return triangles;
}

const std::vector<std::array<int, 2>>& MeshReader::get_BdryEdges_farfield() const
{
    return edgebcfarfield;
}

const std::vector<std::array<int, 2>>& MeshReader::get_BdryEdges_wall() const
{
    return edgebcwall;
}

std::string MeshReader::get_curve_bc_name(int curveEntityTag) const
{
    auto it = curveToPhysicalTags.find(curveEntityTag);
    if (it == curveToPhysicalTags.end())
        return "";

    for (int physicalTag : it->second)
    {
        auto nameIt = physicalNames1D.find(physicalTag);
        if (nameIt != physicalNames1D.end())
            return nameIt->second;
    }

    return "";
}

void MeshReader::read_physical_names(std::ifstream& file)
{
    std::string line;

    std::getline(file, line);
    int nPhysicalNames = std::stoi(line);

    for (int i = 0; i < nPhysicalNames; i++)
    {
        std::getline(file, line);
        std::istringstream iss(line);

        int dim, tag;
        iss >> dim >> tag;

        std::string name;
        iss >> std::ws;
        std::getline(iss, name);

        if (!name.empty() && name.front() == '"') name.erase(0, 1);
        if (!name.empty() && name.back()  == '"') name.pop_back();

        if (dim == 1)
            physicalNames1D[tag] = name;
    }

    std::getline(file, line); // $EndPhysicalNames
}

void MeshReader::read_entities(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);

    std::istringstream header(line);

    int numPoints   = 0;
    int numCurves   = 0;
    int numSurfaces = 0;
    int numVolumes  = 0;

    header >> numPoints >> numCurves >> numSurfaces >> numVolumes;

    // Points
    for (int i = 0; i < numPoints; i++)
    {
        std::getline(file, line);
        // format point:
        // tag x y z numPhysicalTags [physicalTags...] ...
        // on n'en a pas besoin ici
    }

    // Curves
    for (int i = 0; i < numCurves; i++)
    {
        std::getline(file, line);
        std::istringstream iss(line);

        int tag;
        double xmin, ymin, zmin, xmax, ymax, zmax;
        int numPhysicalTags;

        iss >> tag >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax >> numPhysicalTags;

        std::vector<int> phys(numPhysicalTags);
        for (int k = 0; k < numPhysicalTags; k++)
            iss >> phys[k];

        curveToPhysicalTags[tag] = phys;

        // ensuite il reste numBoundingPoints puis leurs tags
        // mais on n'en a pas besoin ici
    }

    // Surfaces
    for (int i = 0; i < numSurfaces; i++)
        std::getline(file, line);

    // Volumes
    for (int i = 0; i < numVolumes; i++)
        std::getline(file, line);

    std::getline(file, line); // $EndEntities
}

