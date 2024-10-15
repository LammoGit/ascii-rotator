#include "graphics.cpp"
#include <span>
#include <vector>
#include <limits>
#include <cstring>
#include <string>
#include <unistd.h>
#include <ncurses.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;

class Screen {
	Mesh mesh;
	Vertice rot;
	Camera camera;
	int width, height, size;
	double *buffZ, spf;
	char *buff;
	public:
		Screen(const Mesh mesh, const Vertice rot, const double fps = 60, const double focus = 30) {
			initscr();
			refresh();
			getmaxyx(stdscr, height, width);

			this->mesh = mesh;
			this->spf = 1.0/fps;
			this->rot = rot * (this->spf * M_PI/180);
			this->camera = Camera(Vertice(0, 0, 50), Vertice(0, 0, -1), width, height, focus);
			this->size = height * width;

			this->buff = new char[size];
			this->buffZ = new double[size];
		}

		~Screen() {
			delete[] buff;
			delete[] buffZ;
			endwin();
		}

		void put(const Vertice &vert) {
			int idx = camera.get_idx(vert);
			if(idx < 0 || idx >= size) return;
			if(buffZ[idx] < camera.distance(vert)) return;
			buff[idx] = '#';
			buffZ[idx] = camera.distance(vert);
		}

		void draw() {
			memset(buff, ' ', size);
			fill(buffZ, buffZ+size, numeric_limits<double>::infinity());
		
			double T, dt = 0.1;
			Vertice k, v;

			for(vector<size_t> &face : mesh.faces) {
				int j;
				for(size_t i=0;i<face.size();i++) {
					j = (i+1) % face.size();
					
					const Vertice &v1 = mesh.vertices.at(face.at(i));
					const Vertice &v2 = mesh.vertices.at(face.at(j));

					k = v2 - v1;
					T = k.length();
					k *= 1/T;

					v = v1;
					for(double t = 0;t < T;t += dt) {
						v += k * dt;
						put(v);
					}
				}
			}
		}

		void run() {
			while(true) {
				clear();
				draw();
				printw(buff);
				refresh();
				mesh.rotate(rot);
				sleep(spf);
			}
		}
};
