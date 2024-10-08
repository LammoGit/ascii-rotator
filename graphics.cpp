#include <span>
#include <format>
#include <math.h>
#include <string>
#include <vector>

template <typename T>
concept Scalar = std::integral<T> || std::floating_point<T>;

double inverse_sqrt(const double x) {
	return 1 / sqrt(x);
}

class Vertice {
	public:
		double x, y, z;
		Vertice(const double x = 0, const double y = 0, const double z = 0) : x{x}, y{y}, z{z} {}
		Vertice(const std::span<Scalar auto> &xyz) {
			if(xyz.size() != 3) { throw format("Span of size {} was given when needed span of size 3.", xyz.size()); }
		}
		inline double length() const { return sqrt(x*x + y*y + z*z); }
		inline double distance(const Vertice &vert) const { return sqrt(pow(x - vert.x, 2)  + pow(y - vert.y, 2)  + pow(z - vert.z, 2)); }
		inline double cosine(const Vertice &vert) const { return (x*vert.x + y*vert.y + z*vert.z)*inverse_sqrt(square_length()*vert.square_length()); }
		inline double project(const Vertice &vert) const { return (x*vert.x + y*vert.y + z*vert.z)*vert.inverse_length(); }

		inline double square_length() const { return x*x + y*y + z*z; }
		inline double inverse_length() const { return inverse_sqrt(x*x + y*y + z*z); }
		inline double square_distance(const Vertice &vert) const { return pow(x - vert.x, 2)  + pow(y - vert.y, 2)  + pow(z - vert.z, 2); }
		inline double inverse_distance(const Vertice &vert) const { return inverse_sqrt(distance(vert)); }

		void rotate(const Vertice &rot) {
			double sx = std::sin(rot.x);
			double sy = std::sin(rot.y);
			double sz = std::sin(rot.z);

			double cx = std::cos(rot.x);	
			double cy = std::cos(rot.y);
			double cz = std::cos(rot.z);
			
			double x0 = x;
			double y0 = y;
			double z0 = z;
			
			x = x0*cy*cz + y0*(sx*sy*cz-cx*sz) + z0*(cx*sy*cz+sx*sz);
			y = x0*cy*sz + y0*(sx*sy*sz+cx*cz) + z0*(cx*sy*sz-sx*cz);
			z = -x0*sy + y0*sx*cy + z0*cx*cy;
		}
		Vertice rotated(const Vertice &rot) const {
			Vertice vert = *this;
			vert.rotate(rot);
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
			double oos = 1/s;
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

		operator std::string() const { return std::format("({}, {}, {})", x, y, z); }
};

class Mesh {
	public:
		std::vector<std::vector<size_t>> faces;
		std::vector<Vertice> vertices;

		Mesh() {}
		Mesh(std::vector<Vertice> vertices, std::vector<std::vector<size_t>> faces) : vertices{vertices}, faces{faces} {}
		
		void rotate(const Vertice &rot) {
			for(Vertice &vert : vertices) {
				vert.rotate(rot);
			}
		}
};

class Camera {
	private:
		Vertice psi, phi;
		Vertice position, direction;
		double half_width, half_height, focus;
		int width;
	public:
		Camera() {}
		Camera(const Vertice position, const Vertice direction, const int width, const int height, const double focus = 30) {
			this->position = position;
			this->width = width;
			this->focus = focus;
			this->half_width = width * 0.5;
			this->half_height = height * 0.5;

			set_direction(direction);
			
		}
		
		const double &get_focus() { return this->focus; }
		const Vertice &get_position() { return this->position; }
		const Vertice &get_direction() { return this->direction; }

		void set_direction(const Vertice direction) {
			if(direction.length() != 1) throw "Length of direction vector should be equal to 1";
			this->direction = direction;
			double sx, sy, cx, cy;
			sy = -direction.x;
			cy = std::sqrt(1 - direction.x * direction.x);
			if(cy) {
				sx = direction.y / cy;
				cx = std::sqrt(1 - sx * sx);
			} else {
				sx = 0;
				cx = 1;
			}
			this->phi = Vertice(cy, -sx*sy, -cx*sy);
			this->psi = Vertice(0, -cx, sx) * 0.5; // I divided psi by 2 because characters in terminal aren't squares
		}

		int get_idx(const Vertice &vert) const {
			Vertice r = vert - position;
			const double nr = r * direction;
			if(nr <= 0) return -1;
			r *= focus / nr;
			const int x = static_cast<int>(r * phi + half_width);
			const int y = width * static_cast<int>(r * psi + half_height);
			return x + y;
		}

		double distance(const Vertice &vert) const {
			return position.square_distance(vert);
		}
};
