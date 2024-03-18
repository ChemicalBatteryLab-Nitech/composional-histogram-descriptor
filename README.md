# Chemical formula derived Histogram descriptor

Histogram descriptor creation program derived from chemical composition formula

Created by Tsubasa KOYAMA (Nagoya Institute of Technology, Nakayama Lab)


## Summary
This is a descriptor creation script that converts various types of information (atomic number, electronegativity, etc.) derived from chemical composition into a general-purpose histogram format. It can be used as a descriptor for machine learning analysis in Materials Infomatics (output is in continuous quantity format).

## detailed explanation
When creating descriptors for a particular system, it is possible to use the corresponding values without adopting a histogram descriptor. However, if the systems are not unified in the target material group, problems such as "defective values are generated in the descriptors" or "descriptors with different meanings are entered in the same column" may occur. If the system is not unified in the material group, problems such as "defective descriptors" or "descriptors with different meanings are entered in the same column" will occur. In machine learning, it is basically difficult to learn when missing values are generated. Also, if descriptors with different meanings are entered in the same column, information that is chemically meaningless for comparison is compared, making it impossible to create a highly accurate prediction model. Histogram descriptors can avoid such problems. Figure 1 shows a model for creating a histogram descriptor using electronegativity (EN) for Li10Zn3Ge4O6 as an example.

![image](https://user-images.githubusercontent.com/106161035/179660726-05805eea-46f3-407f-8a4c-46d5e0ec1325.png)

In Figure 1, electronegativity values are separated at appropriate intervals, and a general vector format descriptor is created by calculating the concentration of the element that falls within each interval. However, since machine learning cannot learn the adjacency of these delimited intervals, an appropriate Gaussian function is applied to smooth the histogram. By representing the multiple information that such compositions contain in vector form, the above problem can be avoided and can be handled for any composition. The elemental properties that can be converted into histogram descriptors with this script are shown in Table 1. Figure 2 shows an example of converting the properties shown in Table 1 into histogram descriptors and combining the vectors. The histogram descriptor is created in this manner.

![image](https://user-images.githubusercontent.com/106161035/179660789-8307643e-cf73-4128-ab5a-0916b501c481.png)
![image](https://user-images.githubusercontent.com/106161035/179660851-be54716f-4e81-47e1-a336-797c11b5581d.png)


## treatment
1. create a csv file with a column listing the chemical formulas (description example: LiCoO2, LiZr2(PO4)3 ) to be converted to descriptors. However, the first row of the column should be the label row.
2. Load the attached .ipynb file with Jupyter notebook or other software, and follow the comments to write the input csv and output csv file names, and then run the program.


## Licensing and citation  (License, Citing)
**License(About License)**　This software is released under the MIT License, see the LICENSE.

**Citation(Citing)**  R. Jalem, M. Nakayama, Y. Noda, T. Le, I. Takeuchi, Y. Tateyama, H. Yamasaki, "A general representation scheme for crystalline solids based on Voronoi-tessellation real feature values and atomic property data", Sci. Technol. Adv. Mater., 19, 231-242 (2018) [DOI: 10.1080/14686996.2018.1439253](https://doi.org/10.1080/14686996.2018.1439253)

## Funding
科研費  19H05815, 20H02436
