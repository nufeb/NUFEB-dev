SHELL = /bin/sh

CC =		g++

# Google Test root dir
GTEST_DIR = ../GTest/googletest

# Google Test preprocessor flags
GTEST = $(GTEST_DIR)/include

# Google test C++ compiler flags
CXXFLAGS = -g -Wall -Wextra

# All Google Test Headers
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
		$(GTEST_DIR)/include/gtest/internal/*.h

# Build gtest.a and gtest_main.a

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)


# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.

all: gtest.a gtest_main.a

gtest-all.o : $(GTEST_SRCS_)
	$(CC) -I${GTEST} -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CC) -I${GTEST} -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

clean:
	rm *.a *.o
