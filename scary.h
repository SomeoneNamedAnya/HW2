#include <bits/stdc++.h>
#include <typeinfo>
#include <cstdlib>
#include <iostream>
#include <string>
#pragma onec;
using namespace std;

class Base {
public:
	virtual void start() = 0;
	virtual void init(vector<vector<char>>&, size_t) = 0;
	virtual void init_m(vector<vector<char>>&, size_t) = 0;
    virtual ~Base(){}
} ;

std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

////////////////////////////////шаблонный класс Fixed ////////////////////////////
template<typename N, int K>
struct Fixed {
    
    constexpr Fixed(int8_t v): v(((N)v) << K) {}
    constexpr Fixed(int16_t v): v(((N)v) << K) {}
    constexpr Fixed(int32_t v): v(((N)v) << K) {}
    constexpr Fixed(int64_t v): v(((N)v) << K) {}
    constexpr Fixed(float f): v(f * (1 << K)) {}
    constexpr Fixed(double f): v(f * (1 << K)) {}
    constexpr Fixed(): v(0) {}

    template<typename N1, int K1>
    constexpr Fixed(const Fixed<N1, K1> & elem) {
        if (K1 > K) {
            v = (elem.v >> (K1 - K));
        } else {
            v = (elem.v << (K - K1));
        }
    }

    template<typename N1, int K1>
    constexpr Fixed<N, K> & operator=(const Fixed<N1, K1> & other) {
        
		if (K1 > K) {
			v = (other.v >> (K1 - K));
		} else {
			v = (other.v << (K - K1));
		}
        
        return *this;
    }
	
    template<typename N1>
    static constexpr Fixed from_raw(N1 x) {
        Fixed<N, K> ret;
        ret.v = x;
        return ret;
    } 

    N v;
	//operator float() const {return  (v / (float) (1 << K));} 
    operator double() const {return  (v / (double) (1 << K));}
    
    //auto operator<=>(const Fixed<N,K>&) const = default;
    template<typename N1>
    bool operator==(const N1 & other)const {
        Fixed<N, K> temp(other);
        return v == temp.v;
    }
    template<typename N1>
    bool operator<(const N1 & other) const{
        Fixed<N, K> temp(other);
        return v < temp.v;
    }
    template<typename N1>
    bool operator>(const N1 & other) const{
        Fixed<N, K> temp(other);
        return (v > temp.v);
    }
    template<typename N1>
    bool operator>=(const N1 & other) const{
        Fixed<N, K> temp(other);
        return (v >= temp.v) ;
    }
    template<typename N1>
    bool operator<=(const N1 & other) const{
        Fixed<N, K> temp(other);
        return (v <= temp.v) ;
    }
    template<typename N1>
    bool operator!=(const N1 & other) const{
        Fixed<N, K> temp(other);
        return (v != temp.v);
    }


};

template<typename N, int K, typename N1>
Fixed<N, K> operator+(Fixed<N, K> a, N1 b) { 
    Fixed<N, K> temp (b);
    return a + temp;
}

template<typename N, int K, typename N1>
Fixed<N, K> operator+(N1 b, Fixed<N, K> a) { 
    Fixed<N, K> temp (b);
    return a + temp;
}

template<typename N, int K, typename N1>
Fixed<N, K> operator-(Fixed<N, K> a, N1 b) { 
    Fixed<N, K> temp (b);
    return a - temp;
}
template<typename N, int K, typename N1>
Fixed<N, K> operator-(N1 b, Fixed<N, K> a) { 
    Fixed<N, K> temp (b);
    return temp - a;
}

template<typename N, int K, typename N1>
Fixed<N, K> operator*(Fixed<N, K> a, N1 b) { 
    Fixed<N, K> temp (b);
    return a * temp;
}
template<typename N, int K, typename N1>
Fixed<N, K> operator*(N1 b, Fixed<N, K> a) { 
    Fixed<N, K> temp (b);
    return a * temp;
}

template<typename N, int K, typename N1>
Fixed<N, K> operator/(Fixed<N, K> a, N1 b) { 
    Fixed<N, K> temp (b);
    return a / temp;
}
template<typename N, int K, typename N1>
Fixed<N, K> operator/(N1 b, Fixed<N, K> a) { 
    Fixed<N, K> temp (b);
    return temp / a;
}


template<typename N, int K, typename N1, int K1>
Fixed<N, K> operator+(Fixed<N, K> a, Fixed<N1, K1> b) {
    if (K1 > K) {
        return Fixed<N, K>::from_raw(a.v + (b.v << (K1 - K)));
    } else {
        return Fixed<N, K>::from_raw(a.v + (b.v >> (K - K1)));
    }
    
}
template<typename N, int K, typename N1, int K1>
Fixed<N, K> operator-(Fixed<N, K> a, Fixed<N1, K1> b) {
    if (K1 > K) {
        return Fixed<N, K>::from_raw(a.v - (b.v << (K1 - K)));
    } else {
        return Fixed<N, K>::from_raw(a.v - (b.v >> (K - K1)));
    }
    
}
template<typename N, int K, typename N1, int K1>
Fixed<N, K> operator*(Fixed<N, K> a, Fixed<N1, K1> b) {
    if (K1 > K) {
        return Fixed<N, K>::from_raw(((int64_t) a.v * (b.v << (K1 - K))) >> K);
    } else {

       return Fixed<N, K>::from_raw(((int64_t) a.v * (b.v >> (K - K1))) >> K);
    }
    
}

template<typename N, int K, typename N1, int K1>
Fixed<N, K> operator/(Fixed<N, K> a, Fixed<N1, K1> b) {
    if (K1 > K) {
        return Fixed<N, K>::from_raw(((int64_t) a.v << K) / (b.v << (K1 - K)));
    } else {
        return Fixed<N, K>::from_raw(((int64_t) a.v << K) / (b.v >> (K - K1)));
    }
    
}


template<typename N, int K, typename N1>
Fixed<N, K>& operator+=(Fixed<N, K>& a, N1 b) { 
    Fixed<N, K> temp (b);
    return a = a + temp;
}
template<typename N, int K, typename N1>
Fixed<N, K>& operator-=(Fixed<N, K>& a, N1 b) { 
    Fixed<N, K> temp (b);
    return a = a - temp;
}
template<typename N, int K, typename N1>
Fixed<N, K>& operator*=(Fixed<N, K>& a, N1 b) { 
    Fixed<N, K> temp (b);
    return a = a * temp;
}
template<typename N, int K, typename N1>
Fixed<N, K>& operator/=(Fixed<N, K>& a, N1 b) { 
    Fixed<N, K> temp (b);
    return a = a / temp;
}

template<typename N, int K>
double& operator/=(double & b, Fixed<N, K> a) { 
    return b = b / (a.v / (double) (1 << K)) ;
}

template<typename N, int K>
float& operator/=(float & b, Fixed<N, K> a) { 
    return b = b / (a.v / (float) (1 << K)) ;
}
template<typename N, int K>
double& operator*=(double & b, Fixed<N, K> a) { 
    return b = b * (a.v / (double) (1 << K)) ;
}

template<typename N, int K>
float& operator*=(float & b, Fixed<N, K> a) { 
    return b = b * (a.v / (float) (1 << K)) ;
}
template<typename N, int K>
double& operator+=(double & b, Fixed<N, K> a) { 
    return b = b + (a.v / (double) (1 << K)) ;
}

template<typename N, int K>
float& operator+=(float & b, Fixed<N, K> a) { 
    return b = b + (a.v / (float) (1 << K)) ;
}
template<typename N, int K>
double& operator-=(double & b, Fixed<N, K> a) { 
    return b = b - (a.v / (double) (1 << K)) ;
}

template<typename N, int K>
float& operator-=(float & b, Fixed<N, K> a) { 
    return b = b - (a.v / (float) (1 << K)) ;
}

template<typename N, int K>
Fixed<N, K> operator-(Fixed<N, K> x) {
    return Fixed<N, K>::from_raw(-x.v);
}

template<typename N, int K>
ostream &operator<<(ostream &out, Fixed<N, K> x) {
    return out << x.v / (double) (1 << K);
}


////////////////////////////////шаблонный класс Simulation///////////////////////////////////////////
template<typename T1, typename T2, typename T3, size_t N, size_t M>
class Simulation : Base{
public:
~Simulation() {}
    void init_m(vector<vector<char>> &te, size_t init_T) {
        cout << "Нужен другой init";
    }
    void init(vector<vector<char>> &te, size_t init_T) {
		Time = init_T;
		cout << N << " " << M << endl;
        rnd.seed(1337);
        rho[' '] = 0.01;
        rho['.'] = 1000;
        g = 0.1;
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < M; j++) {
				field[i][j] = te[i][j];
			}
		}


        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
		
    }

    void start() {
		//cout << "1234567890!!!\n";
        for (size_t i = 0; i < Time; ++i) {
            
            Fixed<int32_t, 16> total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            // Apply forces from p
            memcpy(old_p, p, sizeof(p));
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, (T2)(force / rho[(int) field[x][y]]));
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    cout << field[x] << "\n";
                }
            }
        }
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<Fixed<int32_t, 16>, deltas.size()> tres;
            Fixed<int32_t, 16> sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            Fixed<int32_t, 16> p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, field, p, velocity);
                pp.swap_with(nx, ny, field, p, velocity);
                pp.swap_with(x, y, field, p, velocity);
            }
        }
        return ret;
    }


    Fixed<int32_t, 16> move_prob(int x, int y) {
        Fixed<int32_t, 16> sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    tuple<T3, bool, pair<int, int>> propagate_flow(int x, int y, T2 lim) {
        last_use[x][y] = UT - 1;
        T3 ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = (lim < (T2)(cap - flow) ? lim : (T2)(cap - flow));
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return make_tuple(vp, 1, make_pair(nx, ny));
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return make_tuple(t, prop && end != pair(x, y), end);
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    
    Fixed<int32_t, 16> random01() {
        Fixed<int32_t, 16> a;
        return a.from_raw((rnd() & ((1 << 16) - 1)));
    }

private:
    /*size_t N;
    size_t M; */
    T2 g;
    Fixed<int32_t, 16> rho[256];
    int dirs[N][M]{};
    char field[N][M + 1] {};
    T1 p[N][M]{}, old_p[N][M];
    size_t Time;
    template<typename F>
    struct VectorField {
        array<F, deltas.size()> v[N][M];
        F &add(int x, int y, int dx, int dy, F dv) {
            return get(x, y, dx, dy) += dv;
        }

        F &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    VectorField<T2> velocity{};
    VectorField<T3> velocity_flow{};
    int last_use[N][M]{};
    int UT = 0;

    struct ParticleParams {
        char type;
        T1 cur_p;
        array<T2, deltas.size()> v;

        void swap_with(int x, int y, char field[N][M + 1] , T1 p[N][M], VectorField<T2> &velocity) {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };
    mt19937 rnd;

};

template<typename T1, typename T2, typename T3>
class Simulation<T1, T2, T3, 0, 0> :Base{
public:
public:
    
    void init(vector<vector<char>> &te, size_t init_T) {
		cout << "Нужен init_m\n";
		return;
		
    }
    void init_m(vector<vector<char>> &te, size_t init_T) {
        //cout << "wertgyhuijokl\n";
		Time = init_T;
        N = te.size();
        M = te[0].size();
		//cout << N << " " << M << endl;
        rnd.seed(1337);
        rho[' '] = 0.01;
        rho['.'] = 1000;
        g = 0.1;
 
        velocity.init(N, M);
        velocity_flow.init(N, M);

        field = new char*[N]{};
        dirs = new int*[N]{};
        p = new T1*[N]{};
        old_p = new T1*[N]{};
        last_use = new int*[N]{};
       
		for(size_t i = 0; i < N; i++) {
            field[i] = new char[M + 1]{};
            dirs[i] = new int[M]{};
            last_use[i] = new int[M]{};
            p[i] = new T1[M]{};
            old_p[i] = new T1[M]{};
			for(size_t j = 0; j < M; j++) {
				field[i][j] = te[i][j];
			}
		}


        for (size_t x = 0; x < N; x++) {
            for (size_t y = 0; y < M; y++) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
            
           
        }
		
    }
    void start() {

        for (size_t i = 0; i < Time; ++i) {
            
            Fixed<int32_t, 16> total_delta_p = 0;

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    old_p[i][j] = p[i][j];
                }
            }

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, (T2)(force / rho[(int) field[x][y]]));
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            velocity_flow.clear();
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    cout << field[x] << "\n";
                }
            }
        }
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<Fixed<int32_t, 16>, deltas.size()> tres;
            Fixed<int32_t, 16> sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            Fixed<int32_t, 16> p = random01() * sum;
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, field, p, velocity);
                pp.swap_with(nx, ny, field, p, velocity);
                pp.swap_with(x, y, field, p, velocity);
            }
        }
        return ret;
    }


    Fixed<int32_t, 16> move_prob(int x, int y) {
        Fixed<int32_t, 16> sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    tuple<T3, bool, pair<int, int>> propagate_flow(int x, int y, T2 lim) {
        last_use[x][y] = UT - 1;
        T3 ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = (lim < (T2)(cap - flow) ? lim : (T2)(cap - flow));
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return make_tuple(vp, 1, make_pair(nx, ny));
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return make_tuple(t, prop && end != pair(x, y), end);
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    
    Fixed<int32_t, 16> random01() {
        Fixed<int32_t, 16> a;
        return a.from_raw((rnd() & ((1 << 16) - 1)));
    }

    ~Simulation() {
         for (int i = 0; i < N; i++) {
             //cout << i << endl;
            if (dirs[i] != nullptr) {
                delete[] dirs[i];
            } //cout << i << endl;
            if (field[i] != nullptr) {
                delete[] field[i];
            } //cout << i << endl;
            if (p[i] != nullptr) {
                delete[] p[i];
            } //cout << i << endl;
             if (old_p[i] != nullptr) {
                delete[] old_p[i];
            } //cout << i << endl;
            if (last_use[i] != nullptr) {
                delete[] last_use[i];
            }
            //cout << i << endl;
        }
        if (dirs != nullptr) {
            delete[] dirs;
        }
        if (field != nullptr) {
            delete[] field;
        }
        if (p != nullptr) {
            delete[] p;
        }
        if (old_p != nullptr) {
            delete[] old_p;
        }
        if (last_use != nullptr) {
            delete[] last_use;
        }
            
    }

private:
    size_t N;
    size_t M;
    T2 g;
    Fixed<int32_t, 16> rho[256];
    int** dirs;
    char** field;
    T1 ** p, ** old_p;
    size_t Time;
    template<typename F>
    struct VectorField {
        size_t N, M;
        void init(size_t _N, size_t _M){
            N= _N;
            M = _M;
            v = new array<F, deltas.size()>*[N]{};
            for (int i = 0; i < N; i++) {
                 v[i] = new array<F, deltas.size()>[M]{};
            }
        }
        void clear() {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                 v[i][j] = {};
                }
            }
        }
        ~VectorField() {
            
            for (unsigned i = 0; i < N; i++) {
                if (v[i] == nullptr) continue;
                delete[] v[i];
            }
            if (v != nullptr) {
                delete[] v;
            }
           
        }
        array<F, deltas.size()> ** v;
        F &add(int x, int y, int dx, int dy, F dv) {
            return get(x, y, dx, dy) += dv;
        }

        F &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    VectorField<T2> velocity;
    VectorField<T3> velocity_flow;
    int** last_use;
    int UT = 0;

    struct ParticleParams {
        char type;
        T1 cur_p;
        array<T2, deltas.size()> v;

        void swap_with(int x, int y, char** field , T1 ** p, VectorField<T2> &velocity) {
            swap(field[x][y], type);
            swap(p[x][y], cur_p);
            swap(velocity.v[x][y], v);
        }
    };
    mt19937 rnd;

};

// Макросы для обработки флагов компеляции
#define S(N, M) {N, M}
#define FLOAT make_pair(1, make_pair(0, 0))
#define DOUBLE make_pair(2, make_pair(0, 0))
#define FIXED(N,M) make_pair(3, make_pair(N, M))
#define FAST_FIXED(N,M) make_pair(4, make_pair(N, M))

// size__ = массив всех встречающихся при компиляции размеров поля
#ifdef SIZE 
constexpr pair<int, int> size__[] = {SIZE};
constexpr int sz_size = sizeof(size__) / sizeof(pair<int, int>);
#else
constexpr pair<int, int> size__[0];
constexpr int sz_size = 0;
#endif
// types = массив всех встречающихся при компиляции типов
constexpr pair<int, pair<int, int>> types[] = {TYPES};
constexpr int sz_types = sizeof(types) / sizeof(pair<int, pair<int, int>>);

// map тип <-> реализация шаблона Simulation
map<pair<pair<int, int>, pair<int, int>>, void *> ref_er;
// 4 вложеных цикла
template <int idx1, int idx2, int idx3, int idx4>
 struct it {

    constexpr void doit() {
		 if constexpr (idx4 == sz_size + 1 && sz_size != 0) {
            if (types[sz_types - sz_types] == DOUBLE) {
                ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = new Simulation<double, double, double, size__[0].first, size__[0].second>();
            } else if (types[sz_types - sz_types] == FLOAT) {
                ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = new Simulation<float, float, float, size__[0].first, size__[0].second>();
            } else if (types[sz_types - sz_types].first == 3) {
                if (types[sz_types - sz_types].second.first == 8) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int8_t, types[sz_types - sz_types].second.second>, Fixed<int8_t, types[sz_types - sz_types].second.second>, Fixed<int8_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                } else if (types[sz_types - sz_types].second.first == 16) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int16_t, types[sz_types - sz_types].second.second>, Fixed<int16_t, types[sz_types - sz_types].second.second>, Fixed<int16_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                } else if (types[sz_types - sz_types].second.first == 32) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int32_t, types[sz_types - sz_types].second.second>, Fixed<int32_t, types[sz_types - sz_types].second.second>, Fixed<int32_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                } else if (types[sz_types - sz_types].second.first == 64) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int64_t, types[sz_types - sz_types].second.second>, Fixed<int64_t, types[sz_types - sz_types].second.second>, Fixed<int64_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                }
                
            }   else if (types[sz_types - sz_types].first == 4) {
                if (types[sz_types - sz_types].second.first <= 8) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int_fast8_t, types[sz_types - sz_types].second.second>, Fixed<int_fast8_t, types[sz_types - sz_types].second.second>, Fixed<int_fast8_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                } else if (types[sz_types - sz_types].second.first <= 16) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int_fast16_t, types[sz_types - sz_types].second.second>, Fixed<int_fast16_t, types[sz_types - sz_types].second.second>, Fixed<int_fast16_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                } else if (types[sz_types - sz_types].second.first <= 32) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int_fast32_t, types[sz_types - sz_types].second.second>, Fixed<int_fast32_t, types[sz_types - sz_types].second.second>, Fixed<int_fast32_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                } else if (types[sz_types - sz_types].second.first <= 64) {
                    ref_er[make_pair(make_pair(sz_types - sz_types, sz_types - sz_types),make_pair(sz_types - sz_types, sz_types - sz_types))] = 
                    new Simulation<Fixed<int_fast64_t, types[sz_types - sz_types].second.second>, Fixed<int_fast64_t, types[sz_types - sz_types].second.second>, Fixed<int_fast64_t, types[sz_types - sz_types].second.second>, size__[0].first, size__[0].second>();
                
                }
            }
			if constexpr (idx4 != 0) {
			it<idx1, idx2, idx3, idx4 - 1> x;
			x.doit();
			return;
		} else if constexpr (idx3 != 0) {
			it<idx1, idx2, idx3 - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx2  != 0) {
			it<idx1, idx2 - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx1 != 0) {
			it<idx1 - 1, sz_types - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		}
			
        } else if constexpr (idx4 == sz_size || (sz_size == 0 && idx4 == 1)) {
			
            if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, float, 0, 0>();
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, double, 0, 0>();
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, float, 0, 0>();
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, double, 0, 0>();
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 1){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 2){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 3){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 4){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 1){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 2){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 3){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 4){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, float, 0, 0>();
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, double, 0, 0>();
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, float, 0, 0>();
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, double, 0, 0>();
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 1){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 2){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 3){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 4){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 1){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 2){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 3){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 4){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, float, 0, 0>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, double, 0, 0>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, float, 0, 0>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, double, 0, 0>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, float, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, float, 0, 0>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, double, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, double, 0, 0>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, float, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, float, 0, 0>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, double, 0, 0>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, double, 0, 0>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, 0, 0>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, 0, 0>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, 0, 0>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, 0, 0>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, 0, 0>();
		}
	}
}
}
		if constexpr (sz_size == 0 && idx4 == 1) {
			if constexpr (idx3 != 0) {
				it<idx1, idx2, idx3 - 1, 1> x;
				x.doit();
				return;
			} else if constexpr (idx2 != 0) {
				it<idx1, idx2 - 1, sz_types - 1, 1> x;
				x.doit();
				return;
			} else if constexpr (idx1 != 0) {
				it<idx1 - 1, sz_types - 1, sz_types - 1, 1> x;
				x.doit();
				return;
			}
		} else if constexpr (idx4 != 0) {
			it<idx1, idx2, idx3, idx4 - 1> x;
			x.doit();
			return;
		} else if constexpr (idx3 != 0) {
			it<idx1, idx2, idx3 - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx2  != 0) {
			it<idx1, idx2 - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx1 != 0) {
			it<idx1 - 1, sz_types - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		}
        } else {
           if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, float, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, double, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 1 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, float, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, double, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 2 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 1){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 2){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 3){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 3 && types[idx3].first == 4){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 1){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 2){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 3){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 1 && types[idx2].first == 4 && types[idx3].first == 4){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<float,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, float, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, double, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 1 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 1){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, float, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 2){
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, double, size__[idx4].first, size__[idx4].second>();
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 3){
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 2 && types[idx3].first == 4){
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 1){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 2){
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 3){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 3 && types[idx3].first == 4){
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 1){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 2){
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 3){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 2 && types[idx2].first == 4 && types[idx3].first == 4){
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<double,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 1 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 2 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 3 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 1){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 2){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 3){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 3 && types[idx2].first == 4 && types[idx3].first == 4){
if (types[idx1].second.first == 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, float, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, double, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 1 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,float, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, float, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
} else if (types[idx1].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, double, size__[idx4].first, size__[idx4].second>();
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 2 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
} else if (types[idx1].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,double, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 3 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first == 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 1){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, float, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 2){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	} else if (types[idx2].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, double, size__[idx4].first, size__[idx4].second>();
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 3){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first == 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}
else if (types[idx1].first == 4 && types[idx2].first == 4 && types[idx3].first == 4){
if (types[idx1].second.first <= 8) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast8_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 16) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast16_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 32) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast32_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
} else if (types[idx1].second.first == 64) {
	if (types[idx2].second.first <= 8) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast8_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 16) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast16_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 32) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast32_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	} else if (types[idx2].second.first == 64) {
		if (types[idx3].second.first <= 8) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast8_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 16) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast16_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 32) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast32_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		} else if (types[idx3].second.first == 64) {
			ref_er[{{idx1, idx2}, {idx3, idx4}}] = new Simulation<Fixed<int_fast64_t, types[idx1].second.second>,Fixed<int_fast64_t, types[idx2].second.second>, Fixed<int_fast64_t, types[idx3].second.second>, size__[idx4].first, size__[idx4].second>();
		}
	}
}
}



        if constexpr (idx4 != 0) {
			it<idx1, idx2, idx3, idx4 - 1> x;
			x.doit();
			return;
		} else if constexpr (idx3 != 0) {
			it<idx1, idx2, idx3 - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx2  != 0) {
			it<idx1, idx2 - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		} else if constexpr (idx1 != 0) {
			it<idx1 - 1, sz_types - 1, sz_types - 1, sz_size> x;
			x.doit();
			return;
		}
        }
		
		
  
        
    }

};

template <>
struct it<0, 0, 0, 0> {
    constexpr void doit() {
        return;
    }
};

constexpr void instantiation() { it<sz_types - 1, sz_types - 1, sz_types - 1, sz_size + 1>().doit();}
