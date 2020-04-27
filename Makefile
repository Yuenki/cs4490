all: fps

fps: fps.cpp
	g++ fps.cpp -Wall -lX11 -lGL -lGLU -lm ./libggfonts.a -ofps

clean:
	rm -f fps
	rm -f *.o

