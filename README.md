# Droplet Detective
Droplet Detective aims to bundle the tools and methodologies that are used to perform the so called "slugs test". This is an experimental method that allows for the determination of the tensile strength of concrete and other materials by evaluation the droplet size as the material separates when exiting a nozzle. 

This repository contains:
- a Python script that allows for the continuous logging of load cell measurements
- a Matlab script that allows for the detection of invidual droplets from a log file of a load cell on which droplets have fallen.
- a Matlab script that allows for the estimation of droplet volume from a single image. Based on the assumption that the droplet is axisymmetric.

## References 
- The "Slugs-test" for extrusion-based additive manufacturing: Protocol, analysis and practical limits - https://doi.org/10.1016/j.cemconcomp.2021.104074
- "The Slug Test": Inline Assessment of Yield Stress for Extrusion-Based Additive Manufacturing - https://doi.org/10.1007/978-3-030-49916-7_22
- Automated Visual Inspection of Near Nozzle Droplet Formation for Quality Control of Additive Manufacturing - https://doi.org/10.1007/978-3-031-06116-5_67

## Installation

The current release is tested with Python 3.12 and Matlab R2023b and requires the following libraries:
- To connect to an OPC UA server: opcua-asyncio (py -m pip install asyncua) 

## Coordinates

Gantry robot coordinates for droplet test at TU/e are `X8281 Y4074 Z900`.

## Version numbering
Droplet Detective uses the following [Semantic Versioning](https://semver.org/) scheme: 

```
0.x.x ---> MAJOR version when you make incompatible API changes
x.0.x ---> MINOR version when you add functionality in a backward-compatible manner
x.x.0 ---> PATCH version when you make backward-compatible bug fixes
```

## Contact information

If you have any questions or comments about this project, please open an issue on the repository’s issue page. This can include questions about the content, such as missing information, and the data structure. We encourage you to open an issue instead of sending us emails to help establish an open community. By keeping discussions open, everyone can contribute and see the feedback and questions of others. In addition to this, please see our open science statement below.

## Open science statement

We are committed to the principles of open science to ensure that our work can be reproduced and built upon by others. Our approach includes the following key points:

- Reproducibility: We strive to make our work reproducible by sharing detailed methodologies and data.
- Unique Resources and Data: We have equipment and facilities that are not available at other institutes. We generate data that other institutes cannot produce, and we share this data openly.
- Data and Software Sharing: We share our data and software and encourage others to do the same. If others use our data and software, they must also share their software and data under similar terms.

To support these principles, we license our software under the General Public License version 3 or later (free to use, with attribution, share with source code) and our data and documentation under CC BY-SA (free to use, with attribution, share-alike). By adhering to these principles, we aim to encourage an open and collaborative scientific community. We expect that if you use our resources, you will do the same. 

## License

Copyright (c) 2024 [3D Concrete Printing Research Group at Eindhoven University of Technology](https://www.tue.nl/en/research/research-groups/structural-engineering-and-design/3d-concrete-printing)

Droplet Detective is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License version 3.0 as published by the Free Software Foundation. 

Droplet Detective is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Concrete Candy Tracker; If not, see <http://www.gnu.org/licenses/>.

@license GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.html>
