SRCS := $(wildcard *.cpp)
OBJS := $(addsuffix .o, $(basename $(SRCS)))
GUI_EXE = nufeb_model

all: $(GUI_EXE)

#GeneNet creation Tool Chain
$(GUI_EXE): $(SRCS)
	 g++ nufeb_model.cpp `wx-config-3.0 --cxxflags --libs` -o nufeb_model

.PHONY: clean

clean:
	rm -rf $(OBJS) $(GUI_EXE) *~
