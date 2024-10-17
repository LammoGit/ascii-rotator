#include "graphics.cpp"
#include <span>
#include <vector>
#include <limits>
#include <string>
#include <ncurses.h>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include <unistd.h>
#include <algorithm>

using namespace std;

Vertice get_normal(const Vertice &v1, const Vertice &v2, const Vertice &v3) {
  return (v2-v1).cross((v1-v3).cross(v2-v3)).to_unit();
}

vector<double> intersection(Vertice k1, Vertice b1, Vertice k2, Vertice b2) {
	if(k1 == k2 && b1 == b2) return {0, 1};
	const double det = -k1.x*k2.y+k2.x*k1.y;
	if(det == 0) return {};
	const double oodet = 1/det;
  const double dett = -k2.y*(b2.x-b1.x) + k2.x*(b2.y-b1.y);
	const double dets = k1.x*(b2.y-b1.y) - k1.y*(b2.x-b1.x);
	const double s = dets*oodet;
	if(s >= 0 && s <= 1) return {dett*oodet};
	return {};
}

class Screen {
	Mesh mesh;
	Vertice rot;
  Camera camera;
	LightPoint light;
  int width, height, size;
	double *buffZ, spf;
	char *buff;
	public:
		Screen(const Mesh mesh, const Vertice rot, const LightPoint light, const double fps = 60, const double focus = 30) {
			initscr();
			refresh();
			getmaxyx(stdscr, height, width);

			this->mesh = mesh;
      this->light = light;
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

    char get_char(const Vertice &position) {
      const double brightness = light.brightness(position);
      if(brightness < 0.001) return ' ';
      else if(brightness < 0.01) return '!';
      else if(brightness < 0.2) return '$';
      else if(brightness < 0.4) return '#';
      else return '@';
    }

		void put(const Vertice &vert) {
			int idx = camera.get_idx(vert);
			if(idx < 0 || idx >= size) return;
      double dist = camera.distance(vert);
			if(buffZ[idx] < dist) return;
			buff[idx] = get_char(vert);
			buffZ[idx] = dist;
		}

    void drawline(const Vertice &v1, const Vertice &v2) {
      double T, dt = 0.5;
			Vertice k, v;
			k = v2 - v1;
			T = k.length();
			k *= 1/T;

      v = v1;
			for(double t = 0;t < T;t += dt) {
				v += k * dt;
				put(v);
			}
    }
  
    void flood_lines(const vector<size_t> &face, const Vertice &k, Vertice b, const Vertice &n, const double dt) {
      vector<double> ts, intersections;
      int j;
      do {
        ts.clear();
        for(size_t i=0;i<face.size();i++) {
          j = (i+1) % face.size();

          const Vertice &v1 = mesh.vertices.at(face.at(i));
          const Vertice &v2 = mesh.vertices.at(face.at(j));
          intersections = intersection(k, b, v2-v1, v1);
          ts.insert(ts.end(), intersections.begin(), intersections.end());
        }
        sort(ts.begin(), ts.end());
        for(size_t i=0;i<ts.size();i++) {
          j = (i+1) % ts.size();

          if(ts.at(j) == ts.at(i)) continue;
          const Vertice v1 = k*ts.at(i)+b;
          const Vertice v2 = k*ts.at(j)+b;
          
          drawline(v1, v2);
        }
        b += n*dt;
      } while(ts.size());
    }

    void flood() {
      double dt=0.5;
      Vertice k, b, n;
      for(const vector<size_t> &face : mesh.faces) {
          b = mesh.vertices.at(face.at(0));
          k = mesh.vertices.at(face.at(1)) - b;
          n = get_normal(mesh.vertices.at(face.at(0)), mesh.vertices.at(face.at(1)), mesh.vertices.at(face.at(2)));
          
          flood_lines(face, k, b, n, dt);
          flood_lines(face, k, b, n, -dt);
      }
    }
		
    void draw() {
			memset(buff, ' ', size);
			fill(buffZ, buffZ+size, numeric_limits<double>::infinity());
			flood();
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
