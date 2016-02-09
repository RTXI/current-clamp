PLUGIN_NAME = current_clamp

RTXI_INCLUDES=/usr/local/lib/rtxi_includes

HEADERS = current-clamp.h\
          $(RTXI_INCLUDES)/scatterplot.h\
          $(RTXI_INCLUDES)/incrementalplot.h\
          $(RTXI_INCLUDES)/basicplot.h\
          $(RTXI_INCLUDES)/scrollzoomer.h\
          $(RTXI_INCLUDES)/runningstat.h\
          $(RTXI_INCLUDES)/scrollbar.h\

SOURCES = current-clamp.cpp \
          moc_current-clamp.cpp\
          $(RTXI_INCLUDES)/scatterplot.cpp\
          $(RTXI_INCLUDES)/incrementalplot.cpp\
          $(RTXI_INCLUDES)/basicplot.cpp\
          $(RTXI_INCLUDES)/scrollzoomer.cpp\
          $(RTXI_INCLUDES)/scrollbar.cpp\
          $(RTXI_INCLUDES)/moc_scatterplot.cpp\
          $(RTXI_INCLUDES)/moc_incrementalplot.cpp\
          $(RTXI_INCLUDES)/moc_basicplot.cpp\
          $(RTXI_INCLUDES)/moc_scrollzoomer.cpp\
          $(RTXI_INCLUDES)/moc_scrollbar.cpp\
          $(RTXI_INCLUDES)/runningstat.cpp\
			 
LIBS = -lqwt -lgsl

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
