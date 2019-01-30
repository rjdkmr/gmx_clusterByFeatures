/*
 * This file is part of gmx_clusterByFeatures
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2018  Rajendra Kumar
 *
 * gmx_clusterByFeatures is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * gmx_clusterByFeatures is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with gmx_clusterByFeatures.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>


#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/utility/baseversion.h"

namespace py = pybind11;

int gmx_clusterByFeatures(int argc,char *argv[]);

void wrapped_gmx_clusterByFeatures(std::vector<std::string> argument_vector) {
    char *argv[argument_vector.size()];
    for(size_t n =0; n<argument_vector.size(); n++)
        argv[n] = &argument_vector.at(n)[0];
    
    gmx_run_cmain(argument_vector.size(), argv, &gmx_clusterByFeatures);
}

void wrap_gmx_clusterByFeatures(py::module &m) {
    m.def("gmx_version", &gmx_version);
    m.def("cluster", &wrapped_gmx_clusterByFeatures);
}


PYBIND11_MODULE(gmx_clusterByFeatures, m) 
{
    wrap_gmx_clusterByFeatures(m);
}
