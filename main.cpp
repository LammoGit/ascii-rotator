#include "graphics.cpp"
#include <span>
#include <vector>
#include <limits>
#include <cstring>
#include <string>
#include <unistd.h>
#include <ncurses.h>

using namespace std;

class Screen {
	private:
		vector<Mesh> meshes;
		Vertice rot;
		Camera camera;
		int width, height, size;
		double *buffZ, spf;
		char *buff;

	public:
		Screen(const span<Mesh> meshes, const Vertice rot, const double fps = 60, const double focus = 30) {
			initscr();
			refresh();
			getmaxyx(stdscr, height, width);

			this->meshes.assign(meshes.begin(), meshes.end());
			this->spf = 1.0/fps;
			this->rot = rot * (this->spf * M_PI/180);
			this->camera = Camera(Vertice(0, 0, 50), Vertice(0, 0, -1), width, height, focus);
			this->size = height * width;

			this->buff = new char[size];
			this->buffZ = new double[size];

			run();
		}

		~Screen() {
			delete[] buff;
			delete[] buffZ;
			endwin();
		}

		void draw() {
			memset(buff, ' ', size);
			fill(buffZ, buffZ+size, numeric_limits<double>::infinity());

			for(Mesh &mesh : meshes) {
				draw_mesh(mesh);
			}	
		}

		void put(const Vertice &vert) {
			int idx = camera.get_idx(vert);
			if(idx < 0 || idx >= size) return;
			if(buffZ[idx] < camera.distance(vert)) return;
			buff[idx] = '#';
			buffZ[idx] = camera.distance(vert);
		}

		void draw_mesh(Mesh &mesh) {
			double T, dt = 0.1;
			Vertice k, v;

			for(vector<size_t> &face : mesh.faces) {
				for(size_t &i : face) {
					const Vertice &v1 = mesh.vertices.at(i);
					for(size_t &j : face) {
						if(i == j) continue;
						const Vertice &v2 = mesh.vertices.at(j);
						
						k.x = v2.x - v1.x;
						k.y = v2.y - v1.y;
						k.z = v2.z - v1.z;
						
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
		}

		void rotate() {
			for(Mesh &mesh : meshes) {
				mesh.rotate(rot);
			}
		}

		void run() {
			while(true) {
				clear();
				draw();
				printw(buff);
				refresh();
				rotate();
				sleep(spf);
			}
		}
};

Mesh create_cube(double size, Vertice pos) {
	vector<Vertice> vs;
	vector<vector<size_t>> fs;

	for(int i = -1;i<=1;i+=2)
		for(int j = -1;j<=1;j+=2)
			for(int k = -1;k<=1;k+=2)
				vs.push_back(Vertice(i*size + pos.x, j*size + pos.y, k*size + pos.z));
	fs = { {0, 1, 2, 3}, {0, 1, 4, 5}, {2, 3, 6, 7}, {4, 5, 6, 7}, {0, 2, 4, 6}, {1, 3, 5, 7} };
	return Mesh(vs, fs);
}

Mesh create_pyramid(double size, Vertice pos) {
	vector<vector<size_t>> fs;
	vector<Vertice> vs = { Vertice(size, 0, size), Vertice(size, 0, -size), Vertice(-size, 0, -size), Vertice(-size, 0, size), Vertice(0, size, 0) };
	fs = { {0, 1, 2, 3}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4} };
	return Mesh(vs, fs);
}

int main(int argc, char **argv){
	vector<Mesh> meshes = vector<Mesh>{ create_cube(atof(argv[1]), {0, 0, 0}) };
	Vertice rot = Vertice(atof(argv[2]), atof(argv[3]), atof(argv[4]));
	Screen screen = Screen(meshes, rot, atof(argv[5]), atof(argv[6]));
}
