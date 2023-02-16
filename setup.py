from skbuild import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='baobzi',
      version='0.9.6',
      description='An adaptive fast function approximator based on tree search',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Robert Blackwell',
      author_email='rblackwell@flatironinstitute.org',
      url='https://github.com/flatironinstitute/baobzi',
      packages=['baobzi'],
      package_dir={'baobzi': 'src/python'},
      install_requires=['numpy'],
      cmake_args=['-DBAOBZI_BUILD_TESTS:BOOL=OFF', '-DBAOBZI_BUILD_EXAMPLES:BOOL=OFF', '-DBAOBZI_BUILD_FORTRAN:BOOL=OFF'],
      )
