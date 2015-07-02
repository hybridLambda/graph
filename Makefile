COMPILER = g++
VERSION = ""
FLAGS = -O3 -g -std=c++0x -DVERSION=\"${VERSION}\"
SOURCE = graph.cpp nodeContainer.cpp node.cpp
HEADER = exceptions.hpp graph.hpp nodeContainer.hpp node.hpp

.PHONY: all
all: unittest

unittests_CXXFLAGS = -DUNITTEST -std=c++0x -DVERSION=\"${VERSION}\"
unittests_LDADD = -lcppunit -ldl

test_src = test_graphreader.cpp test_graphbuilder.cpp  ${SOURCE}

unittest: ${test_src} ${HEADER}
	${COMPILER} ${unittests_CXXFLAGS} ${test_src} test_runner.cpp -o unittest ${unittests_LDADD}

.PHONY: clean
clean:
	rm -f unittest

.PHONY: clean_all
clean_all: clean
	rm -f *.o
