# chemhist — Chemical Formula to Histogram Descriptor

A Python package for generating histogram-based compositional descriptors from chemical formulas.

Developed by **Tsubasa Koyama** and **Masanobu Nakayama**  
_Nagoya Institute of Technology, Nakayama Laboratory_

## Overview

`chemhist` converts chemical compositions into histogram-type vector descriptors based on elemental properties (atomic number, electronegativity, etc.).  
These descriptors are designed for use in **Materials Informatics (MI)** and other machine learning analyses, avoiding problems such as missing values or inconsistent feature dimensions among compounds.


## Concept and Motivation

When building descriptors for multiple material systems, direct use of elemental properties may cause problems:
- **Missing values** appear when some elements lack property data.
- **Inconsistent features** occur when descriptors of different meanings occupy the same column.

To solve this, histogram descriptors express elemental properties as continuous distributions.  
Figure 1 illustrates the process using electronegativity (EN) for Li₁₀Zn₃Ge₄O₆.


![image](https://user-images.githubusercontent.com/106161035/179660726-05805eea-46f3-407f-8a4c-46d5e0ec1325.png)


In Figure 1, electronegativity values are separated at appropriate intervals, and a general vector format descriptor is created by calculating the concentration of the element that falls within each interval. However, since machine learning cannot learn the adjacency of these delimited intervals, an appropriate Gaussian function is applied to smooth the histogram. By representing the multiple information that such compositions contain in vector form, the above problem can be avoided and can be handled for any composition. The elemental properties that can be converted into histogram descriptors with this script are shown in Table 1. Figure 2 shows an example of converting the properties shown in Table 1 into histogram descriptors and combining the vectors. The histogram descriptor is created in this manner.

![image](https://user-images.githubusercontent.com/106161035/179660789-8307643e-cf73-4128-ab5a-0916b501c481.png)
![image](https://user-images.githubusercontent.com/106161035/179660851-be54716f-4e81-47e1-a336-797c11b5581d.png)


## Usage

1. Prepare a CSV file containing a column of chemical formulas  
   (e.g., `LiCoO2`, `LiZr2(PO4)3`). The first row must be the header.
2. Open the provided Jupyter notebook (`chemhist_demo.ipynb`) or import the module directly in Python:
   ```python
   from chemhist import get_descriptor
   vec, labels = get_descriptor("Li0.5Mn0.5O2")



## Licensing and citation  (License, Citing)
**License(About License)**　This software is released under the MIT License, see the LICENSE.

**Citation(Citing)**  R. Jalem, M. Nakayama, Y. Noda, T. Le, I. Takeuchi, Y. Tateyama, H. Yamasaki, "A general representation scheme for crystalline solids based on Voronoi-tessellation real feature values and atomic property data", Sci. Technol. Adv. Mater., 19, 231-242 (2018) [DOI: 10.1080/14686996.2018.1439253](https://doi.org/10.1080/14686996.2018.1439253)

## Funding
科研費  19H05815, 20H02436
