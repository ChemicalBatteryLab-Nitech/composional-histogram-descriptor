# chemhist — Chemical Formula to Histogram Descriptor

A Python package for generating histogram-based compositional descriptors from chemical formulas.

Developed by **Tsubasa Koyama** and **Masanobu Nakayama**  
_Nagoya Institute of Technology, Chamical Battery Laboratory_

## Overview

`chemhist` converts chemical compositions into histogram-type vector descriptors based on elemental properties (atomic number, electronegativity, etc.).  
These descriptors are designed for use in **Materials Informatics (MI)** and other machine learning analyses.

## Details

To handle the chemical compositions of materials (mainly inorganic solid compounds) in data science,
it is convenient to represent them as descriptors, which are one-dimensional numerical vectors.
This script converts chemical compositions into histogram descriptors by transforming the elemental properties
that constitute the composition—such as atomic number, electronegativity, and ionic radius—into histograms.
Figure 1 illustrates the generation process of a histogram descriptor using electronegativity (EN)
for the chemical formula Li₁₀Zn₃Ge₄O₆ as an example.


![image](https://user-images.githubusercontent.com/106161035/179660726-05805eea-46f3-407f-8a4c-46d5e0ec1325.png)


The number line shown in Figure 1 represents electronegativity.
By dividing the number line into appropriate intervals and calculating the number of elements contained within each interval
(normalized to a total atomic fraction of 1, i.e., treated as concentration),
a general vector-type descriptor is created.

However, since machine learning models cannot directly learn the adjacency relationships between these discrete intervals,
an appropriate Gaussian function is applied to smooth the histogram,
and the resulting data are output as a vector.

Table 1 lists the elemental properties and their abbreviations that can be converted into histogram descriptors using this script.
Figure 2 illustrates an example in which the properties listed in Table 1 are converted into histogram descriptors
and the resulting vectors are visualized as graphs.<BR>>

**Table 1** Elemental properties used for chemhist descriptors

| Property                     | Abbreviation | Description / Related Quantity     |
|------------------------------|---------------|------------------------------------|
| Electronegativity            | EN            | —                                  |
| Atomic Number                | AN            | —                                  |
| Mendeleev Number             | MN            | —                                  |
| Atomic Weight                | AW            | —                                  |
| Melting Point                | MP            | —                                  |
| Covalent Radius              | CoR           | —                                  |
| Atomic Radius                | AR            | —                                  |
| Ionic Radius                 | IR            | —                                  |
| Crystal Radius               | CrR           | —                                  |
| Group Number                 | PG            | —                                  |
| Period Number                | PN            | —                                  |
| s, p, d, f block elements      | SPDF          | SPDF_0, _1, _2, _3  correspond to concentrations of s-, p-, d-, f- block elements                                  |



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
