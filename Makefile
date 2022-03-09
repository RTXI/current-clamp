PLUGIN_NAME = current_clamp

RTXI_INCLUDES = 

HEADERS = current-clamp.h\

SOURCES = current-clamp.cpp \
          moc_current-clamp.cpp\
			 
LIBS = -lqwt-qt5 -lgsl -lrtplot

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
