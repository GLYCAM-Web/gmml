#include "includes/CentralDataStructure/Readers/Lib/LibraryResidue.hpp"

using lib::LibraryResidue;

LibraryResidue::LibraryResidue(std::stringstream& residueStream, const std::string& name)
{
    this->setName(name);
    this->determineType(name);
    std::string line;
    while(getline(residueStream, line) && line.front() != '!')// Iterate until to the next section that indicates by ! at the beginning of the read line
    { // Process atom section
        //std::cout << "Creating atom with line: " << line << std::endl;
        this->addAtom(std::make_unique<LibraryAtom>(line));
    }
    std::vector<cds::Atom*> createdAtoms = this->getAtoms();
    // Ok atoms are made, now go through the rest of the stream to get their connectivity and coordinates. Not ideal.
    while(getline(residueStream, line) )
    {
        //        if(line.find("boundbox") != std::string::npos)
        //                       double is_box_set;
        //                       std::stringstream ss(line);
        //                       is_box_set << ss ;
        //                       if(is_box_set < 0)              /// If the number written in the first line of the section is negative then boundbox attributes have not been defined in the file
        if(line.find("unit.connect ") != std::string::npos)
        {
            while( getline(residueStream, line) && line.front() != '!')
            {
                /// Process connect section
                std::stringstream ss1(line);
                ss1 >> head_atom_index_;
                getline(residueStream,line);
                std::stringstream ss2(line);
                ss2 >> tail_atom_index_;
            }
        }
        if(line.find("connectivity") != std::string::npos)
        {
            while( getline(residueStream, line) && line.front() != '!')
            {   // Process connectivity section
                std::stringstream ss(line);
                int from, to, t_int;
                ss >> from >> to >> t_int;
                try // yeah this sucks.
                {
                    createdAtoms.at(from-1)->addBond(createdAtoms.at(to-1));
                }
                catch (...)
                {
                    std::stringstream ss;
                    ss << "Tried to access an atom in LibraryResidue constructor that doesn't exist. Positions were " << from << " and " << to << "\n. Number of atoms in residue is " << createdAtoms.size();
                    throw std::runtime_error(ss.str());
                }
            }
        }
        if(line.find("positions") != std::string::npos)
        { // order of atoms in atom vector corresponds to order in this section
            for(auto & atom : createdAtoms)
            {
                getline(residueStream, line);
                std::stringstream ss(line);
                double x, y, z;
                ss >> x >> y >> z;
                cds::Coordinate crd(x, y, z);
                atom->setCoordinate(crd);
            }
        }
    }
    return;
}
