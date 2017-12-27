  Assembly example("example.pdb", PDB);
    example.BuildStructureByDistance(); // Sets the bonding information in AtomNode based on distance
    
    Assembly example1("example.pdb", PDB);
    example1.BuildStructureByDistance(); // Sets the bonding information in AtomNode based on distance

    example.AddSolvent(5, 2, "/home/oliver/Programs/gems/gmml/dat/lib/tip3pbox.off");

    example1.AddSolvent1(5, 2, "/home/oliver/Programs/gems/gmml/dat/lib/tip3pbox.off");
        
    PdbFileSpace::PdbFile *outputPdbFile = example1.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write("old_funciton.pdb");   
    
    PdbFileSpace::PdbFile *outputPdbFile1 = example1.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile1->Write("new_function.pdb");




 solvent_component->BuildAssemblyFromLibraryFile(lib_file);      //Reading the box of water from library file
    solvent_component->SetSourceFile(lib_file);
    solvent_component->BuildStructureByLIBFileInformation();        //Building the structure of the box of water




 Assembly* tip_box = new Assembly();
                tip_box->BuildAssemblyFromLibraryFile(lib_file);
                tip_box->SetSourceFile(lib_file);
                tip_box->BuildStructureByLIBFileInformation();


