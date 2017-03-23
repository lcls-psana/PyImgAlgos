
#------------------------------

import numpy as np
from math import sqrt

#------------------------------

class Quaternion :
    def __init__(self, w=1, x=0, y=0, z=0) :
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def print_q(self) :
        print 'Quaternion w, x, y, z: %.3f  %.3f  %.3f  %.3f' % (self.w, self.x, self.y, self.z)

#------------------------------

class Vector :
    def __init__(self, x, y, z) :
        self.u = self.x = x
        self.v = self.y = y
        self.w = self.z = z

    def print_v(self) :
        print 'Vector: %.3f  %.3f  %.3f' % (self.u, self.v, self.w)

#------------------------------

class Matrix :
    def __init__(self, m00=1, m01=0, n02=0,\
                       m10=0, m11=1, n12=0,\
                       m20=0, m21=0, n22=1) :
        self.m00, self.m01, self.n02 = m00, m01, n02
        self.m10, self.m11, self.n12 = m10, m11, n12
        self.m20, self.m21, self.n22 = m20, m21, n22

    def print_m(self) :
        print 'Matrix:\n   %.3f  %.3f  %.3f\n   %.3f  %.3f  %.3f\n   %.3f  %.3f  %.3f'%\
                 (self.m00, self.m01, self.n02,\
                  self.m10, self.m11, self.n12,\
                  self.m20, self.m21, self.n22)



#------------------------------

def quaternion_modulus(q) :
    """q - quaternion

       If a quaternion represents a pure rotation, its modulus should be unity.
       Returns: the modulus of the given quaternion.
    """
    return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z)

#------------------------------

def normalise_quaternion(q) :
    """q - quaternion

       Rescales the quaternion such that its modulus is unity.

       Returns: the normalised version of q
    """
    mod = quaternion_modulus(q)

    w = q.w / mod
    x = q.x / mod
    y = q.y / mod
    z = q.z / mod

    return Quaternion(w, x, y, z)

#------------------------------

def random_quaternion() :
    """
       Returns: a randomly generated, normalised, quaternion 
    """
    w, x, y, z = 2.0*np.random.random((4,)) - 1.0
    q = Quaternion(w, x, y, z)
    return normalise_quaternion(q)

#------------------------------

def quaternion_valid(q, tol=0.001) :
    """q - quaternion

       Checks if the given quaternion is normalised.

       Returns: 1 if the quaternion is normalised, 0 if not.
    """
    qmod = quaternion_modulus(q)
    if (qmod > 1+tol) or (qmod < 1-tol) : return 0
    return 1

#------------------------------

def quat_rot(v, q) :
    """v - vector (in the form of a "struct rvec") 
       q - quaternion                             
                                                     
       Rotates a vector according to a quaternion.   
                                                     
       Returns: rotated vector vrot            
    """

    t01 = q.w*q.x
    t02 = q.w*q.y
    t03 = q.w*q.z
    t11 = q.x*q.x
    t12 = q.x*q.y
    t13 = q.x*q.z
    t22 = q.y*q.y
    t23 = q.y*q.z
    t33 = q.z*q.z

    u = (1.0 - 2.0 * (t22 + t33)) * v.u\
            + (2.0 * (t12 + t03)) * v.v\
            + (2.0 * (t13 - t02)) * v.w

    v =       (2.0 * (t12 - t03)) * v.u\
      + (1.0 - 2.0 * (t11 + t33)) * v.v\
            + (2.0 * (t01 + t23)) * v.w

    w =       (2.0 * (t02 + t13)) * v.u\
            + (2.0 * (t23 - t01)) * v.v\
      + (1.0 - 2.0 * (t11 + t22)) * v.w

    return Vector(u, v, w)

#------------------------------

def quaternion_from_rotmatrix(m) :
    """m - 3-d rotation matrix, class Matrix 

       Evaluates quaternion from rotation matrix.
       Implemented as in
       https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
       http://www.desy.de/~twhite/crystfel/reference/CrystFEL-Quaternion.html

       Returns: 1 normalised quaternion.
    """
    t, q = Nome, None

    if m.m22 < 0 :
        if m.m00 > m.m11 : 
            t = 1 + m.m00 - m.m11 - m.m22
            q = Quaternion(t, m.m01+m.m10, m.m20+m.m02, m.m12-m.m21)
    
        else : 
            t = 1 - m.m00 + m.m11 - m.m22
            q = Quaternion(m.m01+m.m10, t, m.m12+m.m21, m.m20-m.m02)
    
    else:
        if m.m00 < -m.m11:
            t = 1 - m.m00 - m.m11 + m.m22
            q = Quaternion(m.m20+m.m02, m.m12+m.m21, t, m.m01-m.m10)
    
        else:
            t = 1 + m.m00 + m.m11 + m.m22
            q = Quaternion(m.m12-m.m21, m.m20-m.m02, m.m01-m.m10, t)
    
    #q *= 0.5 / sqrt(t)
    q = normalise_quaternion(q)
    return q

#------------------------------
#------------------------------
#-----------  TEST  -----------
#------------------------------
#------------------------------

def test_quaternion(tname) :

    v1 = Vector(1,0,0)
    v2 = Vector(0,1,0)
    v3 = Vector(0,0,1)
    v1.print_v()
    v2.print_v()
    v3.print_v()

    q1 = Quaternion()
    q1.print_q()

    m1 = Matrix()
    m1.print_m()

#------------------------------

if __name__ == "__main__" :
    import sys; global sys
    tname = sys.argv[1] if len(sys.argv) > 1 else '0'
    print 50*'_', '\nTest %s:' % tname

    test_quaternion(tname)

    sys.exit('End of test %s' % tname)

#------------------------------
