cdef double _h_multiquadric(double r, double sigma):
    return sqrt((1.0/sigma*r)**2 + 1)

cdef double _h_inverse_multiquadric(double r, double sigma):
    return 1.0/sqrt((1.0/sigma*r)**2 + 1)

cdef double _h_thin_plate(double r, double sigma):
    if r == 0.0:
        return 0.0
    else:
        return r**2 * log(r)

cdef double _d_multiquadric(double r, double sigma):
    return r/(sigma**2*sqrt((r/sigma)**2 + 1))

cdef double _d_inverse_multiquadric(double r, double sigma):
    return r/(sigma**2*((r/sigma)**2 + 1)**(3./2))

cdef double _d_thin_plate(double r, double sigma):
    return r * (1 + 2 *log(r))

cdef fitfunction_type setfitfunction(bytes name):

    name = name.lower()
    if name  == "multiquadric":
        funct = _h_multiquadric
    elif name  == "inverse_multiquadric":
        funct =  _h_inverse_multiquadric
    elif name  == "thin_plate":
        funct =  _h_thin_plate
    else:
        raise ValueError("function " + name + " is not supported for fitting")

    return funct

cdef fitfunction_type setfitderivative(bytes name):

    name = name.lower()
    if name  == "multiquadric":
        deriv = _d_multiquadric
    elif name  == "inverse_multiquadric":
        deriv = _d_inverse_multiquadric
    elif name  == "thin_plate":
        deriv = _d_thin_plate
    else:
        raise ValueError("function " + name + " is not supported for fitting")

    return deriv

