a:
	python setup.py build_ext -i
	otool -L Chrono.so
b:
	install_name_tool -change libChronoEngine.dylib @rpath/libChronoEngine.dylib Chrono.so
	otool -L Chrono.so
