CXX=g++
RM=rm -f
CPPFLAGS=-fopenmp -std=c++11
LDFLAGS=-lm
LDLIBS=-fopenmp

SRCS=main.cc test_rig.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(LDFLAGS) -o main_bin $(OBJS) $(LDLIBS)

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS)
	$(RM) *~ .depend

include .depend