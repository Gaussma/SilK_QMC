                      An LA-SiGMA Software Distribution
                                SiLK package
                             Version 1.0 - Jul 2015

The Sign Learning Kink (SiLK) based Quantum Monte Carlo (QMC) method is based
on Feynman's path integral formulation of quantum mechanics, and can reduce the
minus sign problem when calculating energies in atomic and molecular systems.

The code requires as input the one and two electron integrals, which usually
are computed using the NWChem package. Example input files are distributed with
this package. The code also requires an parameter file, specifying run-time
parameters such as input/output directories, or specific code parameters. For
all example inputs a corresponding parameter file is distributed as well.

* Configuration and compilation

  ./src contains the source code and a Makefile. You might need to edit the
  Makefile to fit your system. Needed are suitable C and Fortran compilers,
  and the hdf5, lapack and gsl libraries. The Makefile in this package
  has been tested on a standard Linux installation (Debian Jessie), with
  the relevant libraries installed.

  $ cd src
  $ make

  To clean up the build directory:

  $ make clean

* Running

  To run the code, execute the kink_rwh executable, which needs one command
  line argument: the parameter file. Examples are distributed with this package,
  and are commented inline. They use the file extension '.par'.

  $ ./kink_rwh ./input/H2O_STO3G_c2v/H2O_STO3G_c2v.par

* Example data

  Example data for H2O, F2 and N2 can be found in data/ (also linked as
  src/input), in various basis sets.

Please direct all inquiries about this package to Mark Jarrell
<jarrellphysics@gmail.com>.

This work was supported by the NSF LA-SiGMA EPSCoRprogram  (Cooperative
Agreement No. EPS-1003897) with additional support from the Louisiana Board of
Regents.

LA-SiGMA, the Louisiana Alliance for Simulation-Guided Materials
Applications, is a statewide interdisciplinary collaboration of
material and computer scientists developing computational resources
tuned for the latest generation of accelerator-equipped systems. The
Alliance also develops graduate curricula and is engaged in multiple
outreach activities. Visit us at http://lasigma.loni.org .


The Kink QMC code is an distributed under the terms of the Educational
Community License (ECL) 2.0

Copyright (c) 2012-2013: Louisiana State University
  Xiaoyao Ma <maxiaoyao@gmail.com>
  Frank Löffler <knarf@cct.lsu.edu>
  Randall Hall <randall.hall@dominican.edu>
  Karol Kowalski <karol.kowalski@pnnl.gov>
  Mark Jarrell <jarrellphysics@gmail.com>,
  Juana Moreno <moreno@phys.lsu.edu>

This software and its documentation were developed at Louiana State University. 

Licensed under the Educational Community License, Version 2.0 (the "License"); you may
not use this file except in compliance with the License. You may obtain a copy of the 
License at http://www.osedu.org/licenses/ECL-2.0

Unless required by applicable law or agreed to in writing, software distributed under the 
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, 
either express or implied. See the License for the specific language governing
permissions and limitations under the License.

