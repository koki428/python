if __name__ == '__main__':
    from numpy.distutils.core import setup, Extension

    opt_flags = ['-O3']
    ext = [Extension( name = 'regrid'
                      ,sources = ['regrid.f90'])]
    
    setup(name = 'regrid',
          description = 'Upgrade resolution for R2D2 data',
          author      = 'Hideyuki Hotta',
          ext_modules = ext
          )
