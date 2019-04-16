#!/usr/bin/env python
#
# This file is part of gmx_clusterByFeatures
#
# Author: Rajendra Kumar
#
#
# gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gmx_clusterByFeatures is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.
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

import sys
import binascii

def main():
    pyCode2Hex(sys.argv[1], sys.argv[2])
    
def pyCode2Hex(filename, outfile):
    """ Convert lines in  a file to hex sequences.
    These lines can be read as strings in C++.
    """
    fout = open(outfile, 'w')
    counter = 0
    with open(filename, 'rb') as f:
        for w in f.read():
            if counter == 12:
                fout.write('\n')
                counter = 0
            fout.write("0x{:02x}, ".format(w))
            counter += 1
            
    # Add NULL charecter at the end
    fout.write("0x{:02x}".format(ord('\0')))

if __name__=="__main__":
    main()
