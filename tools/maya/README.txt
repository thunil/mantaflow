Maya 2011 plugins for loading mesh surfaces and density grids
--------------------------------------------------------------

These plugins provide Maya DAG objects which automatically load a 
frame from .bobj.gz or .uni file sequence, corresponding to the current animation frame.
This means that only the current frame has to be held in memory.
The building process has been tested with Maya 2010, 2011 under Windows and Linux,
but should also work with different versions.

Installation
-------------

1. Linux:

> make && make install
This should compile and install the modules in the user's local maya directory.
If maya is installed in a nonstandard directory, you might need to set the MAYA_LOCATION
environment variable or edit buildconfig.linux

2. Windows

TODO

* Download the binaries from the attachement.  Copy bobjFluidObject.mll to a
* directory which is in MAYA_PLUG_IN_PATH (e.g. "My Documents\maya\plug-ins").
* Copy AEbobjFluidObjectTemplate.mel to a directory which is in
* MAYA_SCRIPT_PATH (e.g. "My Documents\maya\scripts").  If you start Maya now,
* bobjFluidObject.mll should be listed in the Plug-in manager in Maya
* ("Windows -> Settings/Preferences -> Plug-in Manager") 

To compile the source code, you'll need zlib. Visual Studio 2008 project files
are included, if you need something else, check the Maya help (under "Developer
Resources -> API Guide -> Setting up your build environment") or use one
of the samples from the Maya devkit as basis.


Usage
------

1. In Maya, you can activate the plugins in Windows->Settings->Plugin Manager.
   If the modules are installed correctly, 'densityloader' and 'bobjloader' should appear.
   Check the load, autoload boxes.
2. To use the objects in a Maya scene, input the following MEL commands:
   > source createBobjLoader; 
   or
   > source createDensityLoader; 
  
3. If you open the attribute editor, you should now be editing the attributes of
   bobjFluidObjectNode1 or gridFluidObjectNode1. 
   Next, you'll need to specify the input file mask. This 
   is the filename of the .bobj.gz-files, containing one printf-like decimal
   placeholder for the frame number (e.g. "C:\data\surface%04d.bobj.gz"). Note
   that the attribute editor script will automatically replace the string "0000"
   with "%04d", so that if you pick the first file of a sequence in the file
   browser, you shouldn't have to make any changes to the file mask.

4. The index offset attribute will be added to the current Maya animation frame
   number to determine the actual file index. The default is -1, because
   .bobj.gz-sequences may be 0-based and Maya frames are 1-based.
   