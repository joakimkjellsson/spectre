import numpy as np
cimport numpy as np

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPE32 = np.int32
ctypedef np.int32_t DTYPE32_t
FTYPE = np.float
ctypedef np.float_t FTYPE_t
FTYPE32 = np.float32
ctypedef np.float32_t FTYPE32_t


cpdef find_bottom_velocity(np.ndarray[np.float32_t, ndim=3] vel, np.ndarray[np.int32_t, ndim=2] kmt, float fill_value=-999.):
   """
   """
   
   cdef int ji,jj,jk,ii,ij,ik
   cdef int imt = vel.shape[2]
   cdef int jmt = vel.shape[1]
   cdef int km = vel.shape[0]
   cdef np.ndarray[np.float32_t, ndim=2] velb = np.ones([jmt,imt], dtype=np.float32)
   
   velb = velb * fill_value
   
   for jj in xrange(0,jmt):
      for ji in xrange(0,imt):
         ik = kmt[jj,ji]
         ik = max(ik,0)
         velb[jj,ji] = vel[ik,jj,ji]
         
   return velb
   