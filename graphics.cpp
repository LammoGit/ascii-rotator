#include <span>
#include <format>
#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <ncurses.h>
#include <stdexcept>

template <typename T>
concept Scalar = std::integral<T> || std::floating_point<T>;

double inverse_sqrt(const double x) {
	return 1 / sqrt(x);
}

class Vertice {
	class VerticeHash {
		public:
			std::size_t operator()(const Vertice &v) const {
				auto h = std::hash<double>();
				return h(v.x) ^ h(v.y) ^ h(v.z);
			}
	};
	inline static std::unordered_map<Vertice, std::vector<double>, VerticeHash> rot_cache;
	public:
		double x, y, z;
		Vertice(const double x = 0, const double y = 0, const double z = 0) : x{x}, y{y}, z{z} {}
		Vertice(const std::span<Scalar auto> &xyz) {
			if(xyz.size() != 3) { throw std::runtime_error(format("Span of size {} was given when needed span of size 3.", xyz.size())); }
		}
		inline double length() const { return sqrt(x*x + y*y + z*z); }
		inline double distance(const Vertice &vert) const { return sqrt(pow(x - vert.x, 2)  + pow(y - vert.y, 2)  + pow(z - vert.z, 2)); }
		inline double cosine(const Vertice &vert) const { return (x*vert.x + y*vert.y + z*vert.z)*inverse_sqrt(square_length()*vert.square_length()); }
		inline double project(const Vertice &vert) const { return (x*vert.x + y*vert.y + z*vert.z)*vert.inverse_length(); }
		inline Vertice cross(const Vertice &vert) const { return Vertice(y*vert.z - z*vert.y, -x*vert.z+z*vert.x, x*vert.y-y*vert.x); }

		inline double square_length() const { return x*x + y*y + z*z; }
		inline double inverse_length() const { return inverse_sqrt(x*x + y*y + z*z); }
		inline double square_distance(const Vertice &vert) const { return pow(x - vert.x, 2)  + pow(y - vert.y, 2)  + pow(z - vert.z, 2); }
		inline double inverse_distance(const Vertice &vert) const { return inverse_sqrt(distance(vert)); }

		static const std::vector<double> &get_rotation(const Vertice &rot) {
			if(rot_cache.contains(rot)) { return rot_cache[rot]; }
			
			const double sx = std::sin(rot.x);
			const double sy = std::sin(rot.y);
			const double sz = std::sin(rot.z);

			const double cx = std::sqrt(1 - sx * sx);
			const double cy = std::sqrt(1 - sy * sy);
			const double cz = std::sqrt(1 - sz * sz);

			rot_cache[rot] = {cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz, cy*sz, sx*sy*sz+cx*cz, cx*sy*sz-sx*cz, -sy, sx*cy, cx*cy};
			return rot_cache[rot];
		}

		void rotate(const Vertice &rot) {
			const std::vector<double> &k = get_rotation(rot);

			const double x0 = x;
			const double y0 = y;
			const double z0 = z;

			x = x0 * k[0] + y0 * k[1] + z0 * k[2];
			y = x0 * k[3] + y0 * k[4] + z0 * k[5];
			z = x0 * k[6] + y0 * k[7] + z0 * k[8];
		}
		void rotate(const std::vector<double> &k) {
			const double x0 = x;
			const double y0 = y;
			const double z0 = z;

			x = x0 * k[0] + y0 * k[1] + z0 * k[2];
			y = x0 * k[3] + y0 * k[4] + z0 * k[5];
			z = x0 * k[6] + y0 * k[7] + z0 * k[8];		
		}
		Vertice rotated(const Vertice &rot) const {
			Vertice vert = *this;
			vert.rotate(rot);
			return vert;
		}
		Vertice rotated(const std::vector<double> &k) const {
			Vertice vert = *this;
			vert.rotate(k);
			return vert;
		}

		Vertice &operator+=(const Scalar auto &s) {
			x += s;
			y += s;
			z += s;
			return *this;
		}
		Vertice &operator-=(const Scalar auto &s) {
			x -= s;
			y -= s;
			z -= s;
			return *this;
		}
		Vertice &operator/=(const Scalar auto &s) {
			const double oos = 1/s;
			x *= oos;
			y *= oos;
			z *= oos;
			return *this;
		}
		Vertice &operator*=(const Scalar auto &s) {
			x *= s;
			y *= s;
			z *= s;
			return *this;
		}

		Vertice &operator+=(const Vertice &vert) {
			x += vert.x;
			y += vert.y;
			z += vert.z;
			return *this;
		}
		Vertice &operator-=(const Vertice &vert) {
			x -= vert.x;
			y -= vert.y;
			z -= vert.z;
			return *this;
		}

		Vertice operator+(const Scalar auto &s) const { return Vertice(x+s, y+s, z+s); }
		Vertice operator-(const Scalar auto &s) const { return Vertice(x-s, y-s, z-s); }
		Vertice operator/(const Scalar auto &s) const { double oos = 1/s; return Vertice(x*oos, y*oos, z*oos); }
		Vertice operator*(const Scalar auto &s) const { return Vertice(x*s, y*s, z*s); }

		Vertice operator+(const Vertice &vert) const { return Vertice(x+vert.x, y+vert.y, z+vert.z); }
		Vertice operator-(const Vertice &vert) const { return Vertice(x-vert.x, y-vert.y, z-vert.z); }
		double operator/(const Vertice &vert) const { return project(vert); }
		double operator*(const Vertice &vert) const { return x*vert.x + y*vert.y + z*vert.z; }

		bool operator==(const Vertice &vert) const { return x==vert.x && y == vert.y && z==vert.z; }
		
		operator std::string() const { return std::format("({}, {}, {})", x, y, z); }
};

class Mesh {
	void check_face(const std::vector<size_t> &face) const{
		if(face.size() == 3) return;
		if(face.size() < 3) throw std::runtime_error("Face must have at least 3 vertices");
		Vertice n = (vertices[face[1]] - vertices[face[0]]).cross(vertices[face[2]] - vertices[face[0]]);
		n /= -(vertices[face[0]] * n);
		for(size_t i=3;i<face.size();i++) {
			if(vertices[face[i]] * n != -1) { 
				throw std::runtime_error("Vertices are on different planes");
			}
		}
	}	
	public:
		std::vector<std::vector<std::size_t>> faces;
		std::vector<Vertice> vertices;

		Mesh() {}
		Mesh(std::vector<Vertice> vertices, std::vector<std::vector<std::size_t>> faces) : vertices{vertices}, faces{faces} {
			if(vertices.size() < 3) throw std::runtime_error("Mesh must have at least 3 vertices");
			for(std::vector<std::size_t> &face : faces) {
				check_face(face);
			}
		}
		
		void rotate(const Vertice &rot) {
			const std::vector<double> &k = Vertice::get_rotation(rot);
			for(Vertice &vert : vertices) {
				vert.rotate(k);
			}
		}
};

class Camera {
	Vertice psi, phi;
	Vertice position, direction;
	double half_width, half_height, focus;
	int width, height;
	public:
		Camera() {}
		Camera(const Vertice position, const Vertice direction, const int width, const int height, const double focus = 30) {
			this->width = width;
			this->height = height;
			this->focus = focus;
			this->half_width = width * 0.5;
			this->half_height = height * 0.5;

			set_direction(direction);
			set_position(position);
		}
		
		const double &get_focus() const { return this->focus; }
		const Vertice &get_position() const { return this->position; }
		const Vertice &get_direction() const { return this->direction; }

		void set_direction(const Vertice direction) {
			this->direction = direction/direction.length();
			double sx, sy, cx, cy;
			sx = direction.y;
			cx = std::sqrt(1-sx*sx);
			if(cx) {
				sy = direction.x/cx;
				cy = -direction.z/cx;
			} else {
				sy = 0;
				cy = 1;
			}
			this->phi = Vertice(cy, 0, sy);
			this->psi = Vertice(-sx*sy, -cx, sx*cy);
		}

		void set_position(const Vertice position) {
			this->position = position;
		}

		int get_idx(const Vertice &vert) const {
			Vertice r = vert - position; // radius vector
			const double nr = r * direction;
			if(nr <= 0) return -1;
			r *= focus / nr;
			const int x = static_cast<int>(r.project(phi) + half_width);
			const int y = static_cast<int>(r.project(psi) * 0.5 + half_height);
			return (x < 0 || y < 0 || x >= width || y >= height) ? -1 : x + y * width;
		}

		double distance(const Vertice &vert) const {
			return position.square_distance(vert);
		}
};
