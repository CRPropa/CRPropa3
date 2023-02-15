# Plugin template

To create own CRPropa modules in C++ the use of plugins is recomended. This plugins can be installed as a seperate module. 
In this folder we provide a template to create such a plugin. 

## general structure of the template
The module class is defined in the header file `myPlugin.h` and the functions in the c++ file `myPlugin.cc`. The new modules are defined in this files. The file `myPlugin.i` defines the SWIG interface for the python usage. 

The python script `testPlugin.py` tests your installation (see below) of the example plugin as presented here.

## adjusting custom module name
To create your own module with a meaningfull modulename all files called `myPlugin.*` and the folder `python/myPlugin` have to be renamed to your plugin name. Also the content of the following files has to be adjusted: 
- `CMakeList.txt`:  The plugin name (see line 4) has to be changed.
- `python/myPlugin/__init__.py`: The folder name has to be changed and change `.myModule` to `.<MyModuleName>`
- `myPlugin.i`: at two positions the header file is included. The name has to be changed. 

# installation of a plugin
For the installation of the plugin you need a running CRPropa version (see [installation documentation](https://crpropa.github.io/CRPropa3/pages/Installation.html)).

First create a build folder and enter it.

    mkdir build && cd build/

Now you can run ccmake to configure your project:

    ccmake ..

At this step you have to set the installation path and the path to your swig interface of the current crpropa installation (if it is not found by cmake).

After configuration (c) and generation (g) you can now build and install your plugin

    make install


## optional testing
Now you can run the python test script. 

    python ../testPlugin.py

