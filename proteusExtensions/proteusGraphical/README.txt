Overview:

The proteusGraphical module is for runtime visualization and rendering
of the solution. It requires Qt and VTK as well as the python wrappers
for Qt and VTK.

Qt/PyQt:

This module requires the Qt 4.5 libraries (not the SDK), PyQt 4.5
(still requires a development snapshot as of 5/8/09), and SIP 4.8
(still requires a development snapshot as of 5/8/09). These are
available from

http://www.qtsoftware.com/downloads
http://www.riverbankcomputing.co.uk/software/sip/download
http://www.riverbankcomputing.co.uk/software/pyqt/download

Install these first following the instructions at their respective
download sites.

VTK:

Go to the external packages directory and install VTK in a directory
called VTK-bin using the source in ParaView3/VTK. You can also install
ParaView3 using the source in ParaView3 but to my knowledge there is
no way to use the VTK libraries that ParaView3 builds for itself. For
convenience a CMakeCache file is provided, but it will require some
modifications for your site.

proteusGraphical: 

python setup.py install

running proteus:

Use the -V vtk switch
