.. _introduction:

Introduction
************

SICOPOLIS (SImulation COde for POLythermal Ice Sheets) is a 3D dynamic/thermodynamic model that simulates the evolution of large ice sheets and ice caps. It was originally created by Greve :cite:`greve_1997a, greve_1997b` in a version for the Greenland ice sheet. Since then, SICOPOLIS has been developed continuously and applied to problems of past, present and future glaciation of Greenland, Antarctica, the entire northern hemisphere, the polar ice caps of the planet Mars and others.

The model employs either hybrid shallow-ice–shelfy-stream dynamics (Bernales et al. :cite:`bernales_etal_2017a, bernales_etal_2017b`) or the shallow-ice approximation for grounded ice, and the shallow-shelf approximation for floating ice (e.g., Greve and Blatter :cite:`greve_blatter_2009`). It is coded in Fortran and uses finite difference discretization on a staggered Arakawa C grid (Arakawa and Lamb :cite:`arakawa_lamb_1977`). A variety of different thermodynamics solvers are available, namely the polythermal two-layer method, two versions of the one-layer enthalpy method, the cold-ice method and the isothermal method (Greve and Blatter :cite:`greve_blatter_2016`).

The coding is based on a low-tech, ease-of-use philosophy. All structures are kept as simple as possible, and advanced coding techniques are only employed where it is deemed appropriate. The use of external libraries is kept at an absolute minimum, which makes the installation very easy and fast.

Resources
=========

* Model website: https://www.sicopolis.net/
* This user manual: https://sicopolis.readthedocs.io/
* GitLab repository: https://gitlab.awi.de/sicopolis/sicopolis/
* SICOPOLIS community @ Zenodo: https://zenodo.org/communities/sicopolis/

Authorship
==========

+-----------------------------------+---------------------------------------------------------+
| Name (affiliation)                | Area of contribution                                    |
+===================================+=========================================================+
| Jorge Bernales (UU)               | Hybrid shallow-ice-shelfy-stream dynamics               |
+-----------------------------------+---------------------------------------------------------+
| Sebastian Beyer (AWI)             | Subglacial hydrology                                    |
+-----------------------------------+---------------------------------------------------------+
| Heinz Blatter (ETHZ)              | Enthalpy method                                         |
+-----------------------------------+---------------------------------------------------------+
| Reinhard Calov (PIK)              | Coupling with climate models, mass balance conservation |
+-----------------------------------+---------------------------------------------------------+
| Liz Curry-Logan (UT Austin)       | Principal developer of SICOPOLIS-AD v1                  |
+-----------------------------------+---------------------------------------------------------+
| Thorben Dunse (HVL)               | Calving, Austfonna                                      |
+-----------------------------------+---------------------------------------------------------+
| Eduardo Flandez (UACh)            | Mocho-Choshuenco ice cap                                |
+-----------------------------------+---------------------------------------------------------+
| Shreyas Sunil Gaikwad (UT Austin) | Principal developer of SICOPOLIS-AD v2                  |
+-----------------------------------+---------------------------------------------------------+
| Ben Galton-Fenzi (UTAS)           | Ice-shelf basal melting                                 |
+-----------------------------------+---------------------------------------------------------+
| Thomas Gölles (UNIS)              | Tracer transport, scripting                             |
+-----------------------------------+---------------------------------------------------------+
| Ralf Greve (HU)                   | Principal developer of SICOPOLIS                        |
+-----------------------------------+---------------------------------------------------------+
| Björn Grieger (ESAC)              | Martian north polar ice cap                             |
+-----------------------------------+---------------------------------------------------------+
| Philipp Hancke (TUD)              | Tibetan Plateau                                         |
+-----------------------------------+---------------------------------------------------------+
| Laurent Hascoet (INRIA)           | Developer of the AD tool Tapenade                       |
+-----------------------------------+---------------------------------------------------------+
| Patrick Heimbach (UT Austin)      | Automatic differentiation (AD)                          |
+-----------------------------------+---------------------------------------------------------+
| Nina Kirchner (SU)                | Tibetan Plateau                                         |
+-----------------------------------+---------------------------------------------------------+
| Thomas Kleiner (AWI)              | Subglacial hydrology                                    |
+-----------------------------------+---------------------------------------------------------+
| Sascha Knell (TUD)                | Elastic lithosphere                                     |
+-----------------------------------+---------------------------------------------------------+
| Anne Le Brocq (UoE)               | Subglacial hydrology                                    |
+-----------------------------------+---------------------------------------------------------+
| Sri Hari Krishna Narayanan (ANL)  | Automatic differentiation (AD)                          |
+-----------------------------------+---------------------------------------------------------+
| Alex Robinson (UCM)               | Coupling with climate models                            |
+-----------------------------------+---------------------------------------------------------+
| Fuyuki Saito (JAMSTEC)            | Code optimization, scripting                            |
+-----------------------------------+---------------------------------------------------------+
| Tatsuru Sato (NIT)                | Ice shelves                                             |
+-----------------------------------+---------------------------------------------------------+
| Marius Schäfer (UACh)             | Northern Patagonian ice field, calving                  |
+-----------------------------------+---------------------------------------------------------+
| Matthias Scheiter (ANU)           | Mocho-Choshuenco ice cap                                |
+-----------------------------------+---------------------------------------------------------+
| Oliver J. Stenzel (MPS)           | Martian north polar ice cap                             |
+-----------------------------------+---------------------------------------------------------+
| Malte Thoma (AWI)                 | Scripting                                               |
+-----------------------------------+---------------------------------------------------------+
| Roland Warner (UTAS)              | Curvilinear coordinates                                 |
+-----------------------------------+---------------------------------------------------------+

| **Affiliations:** 
| ANL: Argonne National Laboratory, ANU: Australian National University, AWI: Alfred Wegener Institute for Polar and Marine Research, ESAC: European Space Astronomy Centre, ETHZ: Swiss Federal Institute of Technology in Zurich, HU: Hokkaido University, HVL: Western Norway University of Applied Sciences, INRIA: National Institute for Research in Digital Science and Technology, JAMSTEC: Japan Agency for Marine-Earth Science and Technology, MPS: Max Planck Institute for Solar System Research, NIT: National Institute of Technology, Ichinoseki College, PIK: Potsdam Institute for Climate Impact Research, SU: Stockholm University, TUD: Darmstadt University of Technology, UACh: Austral University of Chile, UCM: Complutense University of Madrid, UoE: University of Exeter, UNIS: University Centre in Svalbard, UTAS: University of Tasmania, UT Austin: University of Texas at Austin, UU: Utrecht University.

Legal notes
===========

Copyright 2009–2024 SICOPOLIS Authors.

SICOPOLIS is free and open-source software. It can be redistributed and/or modified under the terms of the `GNU General Public License <https://www.gnu.org/licenses/>`__ as published by the Free Software Foundation, either version 3 of the License, or (at the user's option) any later version.

SICOPOLIS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the `GNU General Public License <https://www.gnu.org/licenses/>`__ for more details.

Acknowledgements
================

Thanks a lot for helpful support, suggestions, comments and questions from many colleagues around the world, including those not already listed as SICOPOLIS authors.

Development of SICOPOLIS has been supported by grants/scholarships from

* Alexander von Humboldt Foundation, Germany,
* Federal State of Hesse, Germany,
* German National Academic Foundation (Studienstiftung des deutschen Volkes),
* German Science Foundation (Deutsche Forschungsgemeinschaft DFG),
* Institute of Low Temperature Science, Hokkaido University, Japan,
* Japan Society for the Promotion of Science (JSPS),
* Japanese Ministry of Education, Culture, Sports, Science and Technology (MEXT),
* U.S. Department of Energy, Office of Science,
* U.S. National Science Foundation (NSF).
