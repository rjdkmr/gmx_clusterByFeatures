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

def main():
    options = {'cluster':'Perform clustering using features and extract clustered trajectories',
               'featuresplot' : 'Feature vs Feature plot to check quality of clustering',
               'distmat': 'Calculate avearge, standard-deviation and variance distance-matrix including contact map.',
               'matplot': 'Plot distmat output matrix file',
               'hole'   : 'Calculate channel or cavity radius using hole program',
               'holeplot': 'To calculate average and plot hole output radius file',
               'holefeatures': 'Write radius as a features for clustering',
               'holeclustersplot' : 'To plot or write radius for clusters separately'
              }

    program = 'gmx_clusterByFeatures '
    if len(sys.argv)<=1:
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] not in options:
        print(' ERROR: "{0}" is not an accepted option.\n' .format(sys.argv[1]))
        show_help(options)
        sys.exit(-1)

    if sys.argv[1] == 'cluster':
        from .gmx_clusterByFeatures import cluster
        cluster([program + 'cluster'] + sys.argv[2:])
        
    if sys.argv[1] == 'featuresplot':
        from . import featuresplotcmd
        featuresplotcmd.main()
    
    if sys.argv[1] == 'distmat':
        from .gmx_clusterByFeatures import distmat
        distmat([program + 'distmat'] + sys.argv[2:])
        
    if sys.argv[1] == 'matplot':
        from . import matplot
        matplot.main()
        
    if sys.argv[1] == 'hole':
        from .gmx_clusterByFeatures import hole
        hole([program + 'hole'] + sys.argv[2:])
        
    if sys.argv[1] == 'holeplot':
        from . import holeplot
        holeplot.main()
        
    if sys.argv[1] == 'holefeatures':
        from . import holefeatures
        holefeatures.main()
        
    if sys.argv[1] == 'holeclustersplot':
        from . import holeclustersplot
        holeclustersplot.main()
        
def show_help(options):
    print(' ==============================')
    print(' Usage:')
    print(' gmx_clusterByFeatures <Option>\n')
    print(' ---------------------')
    print(' Use following options:')
    print(' ---------------------\n')

    for tool in options:
        print('\t{0} : {1}\n'.format(tool, options[tool]))

    print(' ==============================')


if __name__=="__main__":
    main()
