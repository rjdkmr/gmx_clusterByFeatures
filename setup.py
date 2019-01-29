#!/usr/bin/env python
#
# This file is part of pyOBabel
#
# Author: Rajendra Kumar
#
# pyOBabel is a python binding to openbabel chemical toolbox (http://openbabel.org).
# Please cite the original publication of the openbabel:
# O'Boyle et al (2011)
# Open Babel: An open chemical toolbox
# Journal of Cheminformatics 2011 3:33
# https://doi.org/10.1186/1758-2946-3-33 .
#
# pyOBabel is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyOBabel is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyOBabel.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#============================================================================

# Always prefer setuptools over distutils
from setuptools import setup, Extension, find_packages
from setuptools.extension import Library
from setuptools.command.build_ext import build_ext, customize_compiler
from setuptools.command.install import install
from distutils.command.build import build
import sys
import setuptools
import os

here = os.path.abspath(os.path.dirname(__file__))

gromacs_flags = None
extensions = None 
        
def check_gromacs_dirs():
    ''' Check for GROMACS directories and flags
    '''
    import pkgconfig
    
    # If gromacs not in pkgconfig, return from here
    if not pkgconfig.exists('libgromacs'):
        return
    
    out = dict()
    
    cppflags_t = pkgconfig.re.split('\s+', pkgconfig.cflags('libgromacs'))
    out['cppflags'] = ''
    out['ldflags'] = []
    out['include'] = []
    out['lib_dirs'] = []
    out['libs'] = []
    
    # Extract include directory and CXXFLAGS
    for flags in cppflags_t:
        if '-I' in flags:
            out['include'].append(flags[2:])
        else:
            out['cppflags'] += flags + ' '
            
    # Extract lib directory and LDFLAGS
    ldflags_t = pkgconfig.re.split('\s+', pkgconfig.libs('libgromacs'))
    for flags in ldflags_t:
        if '-L' in flags:
            out['lib_dirs'].append(flags[2:])
        elif '-l' in flags:
            out['libs'].append(flags[2:])
        else:
            out['ldflags'].append(flags)
 
    
    return out

def extract_gromacs_flags():
    ''' Extract gromacs include, lib and other flags for compilation
    '''
    global  gromacs_flags
    
    # At first check if gromacs is already available in standard path
    gromacs_flags = check_gromacs_dirs()

    # If gromacs is not available at standard path check for GMX_PATH environment variable,
    # add it to pkg-config and then extract GROMACS directories and flags
    if gromacs_flags is None or 'GMX_PATH' in os.environ:
        if 'GMX_PATH' not in os.environ:
            raise LookupError('GMX_PATH environment variable not found...')
        gmx_path = os.environ['GMX_PATH']
        if not os.path.isdir(gmx_path):
            raise LookupError('GROMACS directory {0} not exist...'.format(gmx_path))
        # Check lib name: it could be lib or lib64
        lib_dir = None
        for entry in os.listdir(gmx_path):
            if 'lib' in entry:
                lib_dir = entry
                break
        os.environ['PKG_CONFIG_PATH'] = os.path.join(gmx_path, lib_dir, 'pkgconfig')
        gromacs_flags = check_gromacs_dirs()
    if gromacs_flags is None:
        raise LookupError("gromacs package not found")


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

def get_extensions():
    '''Extensions that are need to be compiled
    '''
    global extensions
    extensions = [ Extension(
        'gmx_clusterByFeatures.gmx_clusterByFeatures',
        [   'src/pywrapper.cpp',
            'src/gmx_clusterbyfeatures.cpp',
            'src/logstream.cpp',
            'src/do_cluster.cpp',
            ],
        include_dirs=[ get_pybind_include(), get_pybind_include(user=True), 
                      'src', ] + gromacs_flags['include'],
        library_dirs=gromacs_flags['lib_dirs'],
        libraries=gromacs_flags['libs'],
        runtime_library_dirs = gromacs_flags['lib_dirs'],
        language='c++',
        extra_link_args= gromacs_flags['ldflags'],
        ),]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'linux':
        c_opts['unix'] += [ '-static-libstdc++'] # Got From https://github.com/pypa/manylinux/issues/118

    def build_extensions(self):
        # Check for -stdlib=libc++ on macos-clang
        if sys.platform == 'darwin':
            # Only in case of clang, so check for this flag
            if has_flag(self.compiler, '-stdlib=libc++'):
                self.c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args += opts
            if sys.platform == 'linux':
                ext.extra_link_args += ['-static-libstdc++']

        # Remove "-Wstrict-prototypes" flags
        customize_compiler(self.compiler)
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass

        build_ext.build_extensions(self)


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


extract_gromacs_flags()
get_extensions()
setup(
    name='gmx_clusterByFeatures',
    version='0.1.0',
    ext_modules=extensions,
    cmdclass={'build_ext': BuildExt},
    install_requires=['pkgconfig>=1.3', 'pybind11>=2.2'],
    entry_points={'console_scripts': [ 'gmx_clusterByFeatures=gmx_clusterByFeatures:main.main',], },
    packages=find_packages(),
    # metadata for upload to pypi
    author = "Rajendra Kumar",	
    author_email = "rjdkmr@gmail.com",
    url = 'https://github.com/rjdkmr/pyOBabel',
    description = "A python binding to openbabel chemical toolbox (http://openbabel.org)",
    long_description = read('README.rst'),
    keywords = ["Molecular Modeling", "Chemoinformatics", "Computational Chemistry", "Computational Drug Design"],
    license = 'GNU General Public License v3 (GPLv3)',
    classifiers = [
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
