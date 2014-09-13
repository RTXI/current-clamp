PLUGIN_NAME = current_clamp

HEADERS = current-clamp.h\
          /usr/local/lib/rtxi_includes/scatterplot.h\
          /usr/local/lib/rtxi_includes/incrementalplot.h\
          /usr/local/lib/rtxi_includes/basicplot.h\
          /usr/local/lib/rtxi_includes/scrollzoomer.h\
			 /usr/local/lib/rtxi_includes/runningstat.h\
			 /usr/local/lib/rtxi_includes/scrollbar.h\

SOURCES = current-clamp.cpp \
			 moc_current-clamp.cpp\
			 /usr/local/lib/rtxi_includes/scatterplot.cpp\
			 /usr/local/lib/rtxi_includes/incrementalplot.cpp\
			 /usr/local/lib/rtxi_includes/basicplot.cpp\
			 /usr/local/lib/rtxi_includes/scrollzoomer.cpp\
			 /usr/local/lib/rtxi_includes/scrollbar.cpp\
			 /usr/local/lib/rtxi_includes/moc_scatterplot.cpp\
			 /usr/local/lib/rtxi_includes/moc_incrementalplot.cpp\
			 /usr/local/lib/rtxi_includes/moc_basicplot.cpp\
			 /usr/local/lib/rtxi_includes/moc_scrollzoomer.cpp\
			 /usr/local/lib/rtxi_includes/moc_scrollbar.cpp\
			 /usr/local/lib/rtxi_includes/runningstat.cpp\
			 
LIBS = -lqwt -lgsl

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
