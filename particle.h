c     particle.h
c     used in up30, darlam
      parameter(npartmax=1000,krelease=1)
      common/specks/nparticl,xrel,yrel,xpart(npartmax),ypart(npartmax)
c     npartmax is maximum number of particles
c     nparticl is current number of particles
c     (xrel, yrel) gives release point (every time step)
c       - provided via (rel_long,rel_lat) in darlam namelist
c     krelease is sigma release level (turned off for krelease=0)
c       - present code assumes particles stay on this given sigma level
