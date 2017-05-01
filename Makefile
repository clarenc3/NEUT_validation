NEUT_ROOT = /vols/build/t2k/cvw09/neut-devel/branches/neut_5.3.6_new_1pi

NEUT_MACHINE = Linux_pc
NEUT_LIB_DIR = $(NEUT_ROOT)/lib/$(NEUT_MACHINE)/
NEUT_CLASS_DIR = $(NEUT_ROOT)/src/neutclass

NEUT_CLASS_LIBS = $(NEUT_CLASS_DIR)/neutpart.so $(NEUT_CLASS_DIR)/neutfsipart.so \
		  $(NEUT_CLASS_DIR)/neutvect.so $(NEUT_CLASS_DIR)/neutfsivert.so \
		  $(NEUT_CLASS_DIR)/neutvtx.so $(NEUT_CLASS_DIR)/neutrootTreeSingleton.so \
		  $(NEUT_CLASS_DIR)/neutnucfsistep.so $(NEUT_CLASS_DIR)/neutnucfsivert.so

NEUT_LIBRARIES = -L$(NEUT_LIB_DIR) -lneutcore -lnuceff -lnuccorrspl -lpartnuck -lskmcsvc -ltauola $(NEUT_CLASS_LIBS)

NEUT_INCLUDES = $(NEUT_CLASS_DIR) -I$(NEUT_ROOT)/include -I$(NEUT_ROOT)/src/pionsmpl

ROOT_INCLUDES = $(shell root-config --incdir)
ROOT_CXX_FLAGS = $(shell root-config --cflags)
ROOT_LIBRARIES = $(shell root-config --glibs --libs)

CXX = g++
CXXFLAGS = -Wall -fPIC $(ROOT_FLAGS) -g
LDFLAGS = -g
SOFLAGS = -shared -Xlinker -zmuldefs

LIBRARIES := $(ROOT_LIBRARIES) $(NEUT_LIBRARIES)
INCLUDES := $(ROOT_INCLUDES) $(NEUT_INCLUDES)



all: 	getKin getKin_nuc checkNEUT

getKin: getKin.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< -I$(NEUT_INCLUDES) -I$(ROOT_INCLUDES) $(LIBRARIES)

getKin_nuc: getKin_nuc.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< -I$(NEUT_INCLUDES) -I$(ROOT_INCLUDES) $(LIBRARIES)

checkNEUT: checkNEUT.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< -I$(NEUT_INCLUDES) -I$(ROOT_INCLUDES) $(LIBRARIES)

clean:
	rm -fv getKin
	rm -fv getKin_nuc
	rm -fv checkNEUT
