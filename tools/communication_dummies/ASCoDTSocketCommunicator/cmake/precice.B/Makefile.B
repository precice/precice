CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS= -Icomponents/java -Isrc -I/usr/lib/jvm/java-6-openjdk/include -I/usr/lib/jvm/java-6-openjdk/include/linux -fPIC
LDFLAGS=
LDLIBS=

SRCS= components/java/__PATH_FULL_QUALIFIED_NAME__AbstractImplementation.cpp src/__PATH_FULL_QUALIFIED_NAME__Implementation.cpp components/java/__PATH_FULL_QUALIFIED_NAME__JNIProxy.cpp 
SRCS += components/java/precice/InitializerNative2JavaPlainPort.cpp components/java/precice/InitializerNative2NativePlainPort.cpp components/java/precice/InitializerNativeDispatcher.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: libB.so

libB.so: $(OBJS)
	$(CXX) -shared $(OBJS) -o native/libB.so

clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) libB.so
