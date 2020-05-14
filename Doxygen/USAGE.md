# Building the Doxygen documentation

**Important**

1.  The contents of this directory are meant to be used from the directory above. 
    * All the references assume that you are in the directory above.
    * The documentation that is created is placed in this directory.
2.  The Doxyfile in this directory contains variables that need to be assigned.
    * These are assigned in the script:  Build\_Doxygen\_Docs.bash
    * If you don't want to use that script, you can export them 
      yourself on the command line or in your own script.

**Usage**

The preferred way to use the documentation is:


```
bash ./Doxygen/Build_Doxygen_Docs.bash
```

This will generate documentation in the Doxygen/html directory.

To view the files, point a browser window to:

```
/path/to/code/Doxygen/html/index.html
```

