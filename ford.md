project: SICOPOLIS
summary: SICOPOLIS Developer Manual
author: SICOPOLIS Authors
src_dir: ./src
preprocess: false
docmark: !
docmark_alt: *
predocmark: >
predocmark_alt: |
display: public
         protected
         private
incl_src: true
lower: false
max_frontpage_items: 20
search: true
sort: src
source: false
version: v24
year: 2024
graph: true
output_dir: ./docs/ford
force: true

SICOPOLIS (SImulation COde for POLythermal Ice Sheets) is a 3D dynamic/thermodynamic model that simulates the evolution of large ice sheets and ice caps. It was originally created by Greve (1997a,b) in a version for the Greenland ice sheet. Since then, SICOPOLIS has been developed continuously and applied to problems of past, present and future glaciation of Greenland, Antarctica, the entire northern hemisphere, the polar ice caps of the planet Mars and others.

The model employs either hybrid shallow-ice-shelfy-stream dynamics (Bernales et al. 2017a,b) or the shallow-ice approximation for grounded ice, and the shallow-shelf approximation for floating ice (e.g., Greve and Blatter 2009). It is coded in Fortran and uses finite difference discretization on a staggered Arakawa C grid. A variety of different thermodynamics solvers are available, namely the polythermal two-layer method, two versions of the one-layer enthalpy method, the cold-ice method and the isothermal method (Greve and Blatter 2016).

The coding is based on a low-tech, ease-of-use philosophy. All structures are kept as simple as possible, and advanced coding techniques are only employed where it is deemed appropriate. The use of external libraries is kept at an absolute minimum, which makes the installation very easy and fast.

##### References

See <https://sicopolis.readthedocs.io/en/latest/references.html>.

##### Resources

- Model website: <https://www.sicopolis.net/>  
- User manual: <https://sicopolis.readthedocs.io/>  
- GitLab repository: <https://gitlab.awi.de/sicopolis/sicopolis/>  
- SICOPOLIS community @ Zenodo: <https://zenodo.org/communities/sicopolis/>

##### Copyright

Copyright 2009-2024 [SICOPOLIS Authors](https://sicopolis.readthedocs.io/en/latest/introduction.html#authorship).

##### License

SICOPOLIS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SICOPOLIS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SICOPOLIS. If not, see <https://www.gnu.org/licenses/>.

@note
This developer manual is generated automatically using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).  
Installing FORD: `pip install ford`.  
Creating/updating the developer manual: `ford ford.md`.  
Output is in HTML -> `docs/ford/index.html`.
