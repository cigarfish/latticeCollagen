
QT += opengl

macx:LIBS += -framework GLUT
linux-g++:LIBS += -lglut -lGLU
linux-g++-64:LIBS += -lglut -lGLU
win32:LIBS += -lglut32
