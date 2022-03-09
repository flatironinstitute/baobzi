from skbuild import setup

setup(name='baobzi',
      version='0.9.3',
      description='An adaptive fast function approximator based on tree search',
      author='Robert Blackwell',
      author_email='rblackwell@flatironinstitute.org',
      url='https://github.com/flatironinstitute/baobzi',
      packages=['baobzi'],
      package_dir={'baobzi': 'src/python'},
      install_requires=['numpy'],
      cmake_args=['-DBAOBZI_BUILD_TESTS:BOOL=OFF', '-DBAOBZI_BUILD_EXAMPLES:BOOL=OFF', '-DBAOBZI_BUILD_FORTRAN:BOOL=OFF'],
      )
