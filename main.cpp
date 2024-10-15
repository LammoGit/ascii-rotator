#include "screen.cpp"
#include <vector>
#include <cstdlib>

using namespace std;

Mesh create_cube(double size, Vertice pos) {
	vector<Vertice> vs;
	vector<vector<size_t>> fs;

	for(int i=-1;i<=1;i+=2) {
		for(int j=-1;j<=1;j+=2) {
			for(int k=-1;k<=1;k+=2) {
				vs.push_back(Vertice(i*size + pos.x, j*size + pos.y, k*size + pos.z));
			}
		}
	}
	fs = { {0, 1, 5, 4}, {2, 3, 7, 6}, {1, 3, 7, 5}, {0, 2, 6, 4}, {1, 3, 2, 0}, {4, 5, 7, 6} };
	return Mesh(vs, fs);
}

Mesh create_pyramid(double size, Vertice pos) {
	vector<vector<size_t>> fs;
	vector<Vertice> vs = { Vertice(size, 0, size), Vertice(size, 0, -size), Vertice(-size, 0, -size), Vertice(-size, 0, size), Vertice(0, size, 0) };
	fs = { {0, 1, 2, 3}, {0, 3, 4}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4} };
	return Mesh(vs, fs);
}

int main(int argc, char **argv){
	Mesh mesh = create_pyramid(atof(argv[1]), {0, 0, 0});
	Vertice rot = Vertice(atof(argv[2]), atof(argv[3]), atof(argv[4]));
	Screen screen = Screen(mesh, rot, atof(argv[5]), atof(argv[6]));
	screen.run();
}
