# Plugin template

To create your own CRPropa modules in C++ we recommend using plugins. Plugins are small programs that can be installed as a separate modules. 
Here, we provide a template to create such a plugin. 

## General structure of the template
The module class is defined in the header file `myPlugin.h` and the functions in the C++ file `myPlugin.cc`. The file `myPlugin.i` defines the SWIG interface for python usage. 

The python script `testPlugin.py` tests the installation (see below) of the example plugin as presented here.

## Adjusting custom module name
To create your own module with a meaningful module name, all files called `myPlugin.*` and the folder `python/myPlugin` have to be renamed to your plugin name. Also the content of the following files has to be adjusted: 
- `CMakeList.txt`:  The plugin name (see line 4) has to be changed.
- `python/myPlugin/__init__.py`: The directory name and the content of the init-file have to be changed: `.myModule` to `.<MyModuleName>`
- `myPlugin.i`: at two positions the header file is listed. The lines have to be adapted accordingly. 

# Installation of a plugin
For the installation of the plugin you need a running CRPropa version (see [installation documentation](https://crpropa.github.io/CRPropa3/pages/Installation.html)).
This is done analogously to the installation of CRPropa. We recommend to activate the same virtual python environment that you use to run CRPropa.

First create a build folder within the plugin's directory and move there.

    mkdir build && cd build/

Now you can run cmake to configure your project:

    ccmake ..

At this step you have to set the installation path and the path to your swig interface of the current CRPropa installation (if it is not found by cmake).

After configuration (press c) and generation (press g) you can now build and install your plugin

    make install


## optional testing
Now you can run the python test script. 

    python ../testPlugin.py