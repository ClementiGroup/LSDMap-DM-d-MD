ctypedef double (*fitfunction_type)(double, double)

cdef double _h_multiquadric(double r, double sigma)
cdef double _h_inverse_multiquadric(double r, double sigma)
cdef double _h_thin_plate(double r, double sigma)
cdef double _h_gaussian(double r, double sigma)

cdef double _d_multiquadric(double r, double sigma)
cdef double _d_inverse_multiquadric(double r, double sigma)
cdef double _d_thin_plate(double r, double sigma)
cdef double _d_gaussian(double r, double sigma)

cdef fitfunction_type setfitfunction(bytes name)
cdef fitfunction_type setfitderivative(bytes name)

cdef extern from "math.h":
    double exp(double x)
    double sqrt(double x)
    double log(double x)
