/*Field is a container for mesh node data division by volume*/
#ifndef _FIELD_H
#define _FIELD_H

#include <ostream>
#include "Arena.h"

template <typename T>
struct vec3 {
	vec3 (const T u, const T v, const T w) : d{u,v,w} {}
	vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
	vec3 (): d{0,0,0} {}
	T& operator[](int i) {return d[i];}
	T operator()(int i) const {return d[i];}
	vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
	vec3<T>& operator+=(vec3<T> o) {d[0]+=o[0];d[1]+=o[1];d[2]+=o[2];return(*this);}
	vec3<T>& operator-=(vec3<T> o) {d[0]-=o[0];d[1]-=o[1];d[2]-=o[2];return(*this);}

protected:
	T d[3];
};

//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2));	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2));	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2));	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2));	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}
template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;


template <typename T>
class Field_
{
public:
	
	/*constructor*/
	Field_(ArenaAllocator &arena, int ni, int nj, int nk) :
	ni{ni}, nj{nj}, nk{nk}
	{
		//allocate memory for a 3D array
		data = arena.alloc<T>(ni*nj*nk);	

		clear();
	}
	
	/*Fast Indexing Helper*/
	inline int idx(int i, int j, int k) const {
		return i* nj * nk + j * nk + k;
	}

	T operator() (int i, int j, int k) const {return data[idx(i,j,k)];}
	T& operator() (int i, int j, int k) {return data[idx(i,j,k)];}

	/*sets all values to some scalar*/
	void operator =(double s) {
		int max = ni*nj*nk;
		for (int x = 0; x< max; x++) data[x] = s;
	  }

	/*performs element by element division by another field*/
	void operator /= (const Field_ &other) {
		int max = ni*nj*nk;
		for (int x = 0; x< max; x++) {
			if (other.data[x] != 0) data[x] /= other.data[x];
			else data[x] = 0;
		}
	}

	/*increments values by data from another field*/
	Field_& operator += (const Field_ &other) {
		int max = ni*nj*nk;
		for (int x = 0; x< max; x++) data[x] += other.data[x];
		return *this;
	}

	/*performs element by element division by another field*/
	Field_& operator *= (double s) {
		int max = ni * nj * nk;
        for (int x = 0; x < max; x++) data[x] *= s;
        return *this;
	}


	/*sets all data to zero*/
	void clear() {(*this)=0;}
	
	/* 4. Scatter / Gather */
    void scatter(double3 lc, T value) {
        int i = (int)lc[0]; double di = lc[0]-i;
        int j = (int)lc[1]; double dj = lc[1]-j;
        int k = (int)lc[2]; double dk = lc[2]-k;
        
        (*this)(i,   j,   k)   += (T)value*(1-di)*(1-dj)*(1-dk);
        (*this)(i+1, j,   k)   += (T)value*(di)*(1-dj)*(1-dk);
        (*this)(i+1, j+1, k)   += (T)value*(di)*(dj)*(1-dk);
        (*this)(i,   j+1, k)   += (T)value*(1-di)*(dj)*(1-dk);
        (*this)(i,   j,   k+1) += (T)value*(1-di)*(1-dj)*(dk);
        (*this)(i+1, j,   k+1) += (T)value*(di)*(1-dj)*(dk);
        (*this)(i+1, j+1, k+1) += (T)value*(di)*(dj)*(dk);
        (*this)(i,   j+1, k+1) += (T)value*(1-di)*(dj)*(dk);     
    }

    T gather(double3 lc) {
        int i = (int)lc[0]; double di = lc[0]-i;
        int j = (int)lc[1]; double dj = lc[1]-j;
        int k = (int)lc[2]; double dk = lc[2]-k;
        
        return  (*this)(i,   j,   k)   * (1-di)*(1-dj)*(1-dk)+
                (*this)(i+1, j,   k)   * (di)*(1-dj)*(1-dk)+
                (*this)(i+1, j+1, k)   * (di)*(dj)*(1-dk)+
                (*this)(i,   j+1, k)   * (1-di)*(dj)*(1-dk)+
                (*this)(i,   j,   k+1) * (1-di)*(1-dj)*(dk)+
                (*this)(i+1, j,   k+1) * (di)*(1-dj)*(dk)+
                (*this)(i+1, j+1, k+1) * (di)*(dj)*(dk)+
                (*this)(i,   j+1, k+1) * (1-di)*(dj)*(dk);
    }

    // Templated friend declaration for output
    template<typename S>
    friend std::ostream& operator<<(std::ostream &out, Field_<S> &f);

protected:
    int ni, nj, nk;
    T* data; 
};

/* Output Stream Operator */
template<typename T>
std::ostream& operator<<(std::ostream &out, Field_<T> &f)
{
    for (int k=0; k<f.nk; k++, out<<"\n")
        for (int j=0; j<f.nj; j++)
            for (int i=0; i<f.ni; i++) 
                out << f(i,j,k) << " ";
    return out;
}

using Field = Field_<double>;
using Field3 = Field_<double3>;

#endif