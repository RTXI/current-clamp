PLUGIN_NAME = CClamp4

HEADERS = CClamp4.h\
    include/runningstat.h\
    include/basicplot.h\
    include/scatterplot.h\
    include/incrementalplot.h\
    include/scrollbar.h\
    include/scrollzoomer.h\
    include/RTXIprintfilter.h

SOURCES = CClamp4.cpp \
    moc_CClamp4.cpp\
    include/runningstat.cpp\
    include/basicplot.cpp\
    include/scatterplot.cpp\
    include/incrementalplot.cpp\
    include/scrollbar.cpp\
    include/scrollzoomer.cpp\
    include/moc_runningstat.cpp\
    include/moc_scatterplot.cpp\
    include/moc_incrementalplot.cpp\
    include/moc_scrollbar.cpp\
    include/moc_scrollzoomer.cpp\
    include/moc_basicplot.cpp

LIBS = -lqwt -lgsl

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
