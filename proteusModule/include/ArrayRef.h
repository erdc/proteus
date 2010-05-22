#ifndef ARRAYREF_H
#define ARRAYREF_H

template <class T>
class ArrayRef1D
{
	T *values;

public:
	inline ArrayRef1D(T* vals) : values(vals) {};
	inline T &operator[](int i) { return values[i]; };
};

template <class T>
class ArrayRef2D
{
	T *values;
	int dimension1;
	int dimension2;

public:
	inline ArrayRef2D(T* vals, int dim1, int dim2) : values(vals), dimension2(dim2) {};
	inline ArrayRef2D(T* vals, int dim2) : values(vals), dimension2(dim2) {};
	inline T* operator[] (int i) { return values + i * dimension2; }
};

template <class T>
class ArrayRef3D
{
	T *values;
	int dimension1;
	int dimension2;
	int dimension3;
	int mult;

public:
	ArrayRef3D(T* vals, int dim1, int dim2, int dim3) : values(vals), dimension2(dim2), dimension3(dim3), mult(dim2 * dim3) {};
	inline ArrayRef2D<T> operator[] (int i) { return ArrayRef2D<T>(values + i * mult, dimension3); } ;
};

template <class T>
class ArrayRef4D
{
	T *values;
	int dimension1;
	int dimension2;
	int dimension3;
	int dimension4;
	int mult1;

public:
	ArrayRef4D(T* vals, int dim1, int dim2, int dim3, int dim4) : 
	  values(vals), dimension2(dim2), dimension3(dim3), dimension4(dim4),mult1(dim2 * dim3 * dim4) {};
	inline ArrayRef3D<T> operator[] (int i) { return ArrayRef3D<T>(values + i * mult1 , dimension2, dimension3, dimension4); } ;
};

template <class T>
class ArrayRef5D
{
	T *values;
	int dimension1;
	int dimension2;
	int dimension3;
	int dimension4;
	int dimension5;
	int mult1;

public:
	ArrayRef5D(T* vals, int dim1, int dim2, int dim3, int dim4, int dim5) : 
	  values(vals), dimension2(dim2), dimension3(dim3), dimension4(dim4), dimension5(dim5),mult1(dim2 * dim3 * dim4* dim5) {};
	inline ArrayRef4D<T> operator[] (int i) { return ArrayRef4D<T>(values + i * mult1 , dimension2, dimension3, dimension4, dimension5); } ;
};

template <class T>
class ArrayRef6D
{
	T *values;
	int dimension1;
	int dimension2;
	int dimension3;
	int dimension4;
	int dimension5;
	int dimension6;
	int mult1;

public:
	ArrayRef6D(T* vals, int dim1, int dim2, int dim3, int dim4, int dim5, int dim6) : 
	  values(vals), dimension2(dim2), dimension3(dim3), dimension4(dim4), dimension5(dim5), dimension6(dim6), mult1(dim2 * dim3 * dim4 * dim5 * dim6) {};
	inline ArrayRef5D<T> operator[] (int i) { return ArrayRef5D<T>(values + i * mult1 , dimension2, dimension3, dimension4, dimension5, dimension6); } ;
};

#endif
