# chemhist — Chemical Formula to Histogram Descriptor

A Python package for generating histogram-based compositional descriptors from chemical formulas (chemhist descriptors).

Developed by **Tsubasa Koyama** and **Masanobu Nakayama**  
_Nagoya Institute of Technology, Chamical Battery Laboratory_

## 1. Overview

`chemhist` converts chemical compositions into histogram-type vector descriptors based on elemental properties (atomic number, electronegativity, etc.).  
These descriptors are designed for use in machine learning analyses for materials (Materials Informatics).

## 2. Description

### 2.1 Standard Chemhist descriptors
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

![image](https://user-images.githubusercontent.com/106161035/179660851-be54716f-4e81-47e1-a336-797c11b5581d.png)

### 2.2 Algebric Descriptors

### 2.3 Matrix Descriptors


## 3. Installation

`chemhist` can be installed locally from source using `pip`.  
Make sure you are in the directory that contains the `pyproject.toml` file.

```bash
cd path/to/chemhist_project
pip install .
```

If you encounter any build errors (for example, *access denied* or *failed to build wheel*),  
clean up previous build directories and try again:

```bash
# Windows PowerShell
Remove-Item -Recurse -Force build, dist, chemhist.egg-info
```

To verify that the installation was successful:

```python
import chemhist
print(chemhist.__file__)
```

If the package is installed correctly, the path to  
`site-packages/chemhist/__init__.py` will be displayed.


## 4 Usage


### 4.1. Import and generate descriptors
The main function is `get_descriptor()`, which converts a chemical formula into a histogram-based descriptor vector.

```python
from chemhist import get_descriptor

# Example: create histogram descriptor for Li0.5Mn0.5O2
vec, labels = get_descriptor("Li0.5Mn0.5O2")

print("Number of features:", len(vec))
print("First 10 features:", vec[:10])
```

This returns:
- `vec` : a NumPy 1D array containing the descriptor values (broadened histogram)
- `labels` : a list of feature names corresponding to each element of `vec`

---

### 4.2. Generate descriptors for multiple compositions

You can loop through multiple chemical compositions to build a dataset.

```python
from chemhist import get_descriptor
import pandas as pd

formulas = ["LiCoO2", "LiMnO2", "Na3PS4"]

records = []
for f in formulas:
    vec, labels = get_descriptor(f)
    records.append([f] + list(vec))

df = pd.DataFrame(records, columns=["formula"] + labels)
df.to_csv("chemhist_output.csv", index=False)

print("Descriptors saved to chemhist_output.csv")
```

This will create a CSV file like:
```
formula, EN_1, EN_2, EN_3, ..., PG_10
LiCoO2, 0.12, 0.08, 0.00, ..., 0.03
LiMnO2, 0.10, 0.05, 0.02, ..., 0.01
...
```

---

### 4.3. Plot the histogram descriptor

You can visualize the resulting histogram descriptors as a bar chart.

```python
import matplotlib.pyplot as plt
from chemhist import get_descriptor

vec, labels = get_descriptor("Li0.5Mn0.5O2")

plt.figure(figsize=(12,4))
plt.bar(range(len(vec)), vec)
plt.xlabel("Descriptor index")
plt.ylabel("Value")
plt.title("Histogram Descriptor for Li0.5Mn0.5O2")
plt.tight_layout()
plt.show()
```

If you wish to group the bars by property type (e.g., EN, PG, AR...),  
you can use the following snippet:

```python
props = [l.split("_")[0] for l in labels]
unique_props = []
for p in props:
    if p not in unique_props:
        unique_props.append(p)

plt.figure(figsize=(12,4))
for p in unique_props:
    idx = [i for i, l in enumerate(labels) if l.startswith(p)]
    plt.bar(idx, vec[idx], label=p)

plt.legend(ncol=4)
plt.xlabel("Descriptor Index")
plt.ylabel("Value")
plt.title("Histogram Descriptors Grouped by Property")
plt.tight_layout()
plt.show()
```

---

### 4.4. CLI execution (optional)
After installation, you can also run chemhist directly from the command line:

```bash
python -m chemhist LiCoO2 LiMnO2 --out descriptors.csv
```

This command will:
- Convert the listed chemical formulas into histogram descriptors
- Save the results into a CSV file (`descriptors.csv` by default)
- Include the feature names in the header


---

### 4.6 Notes
- The broadening is applied by a Gaussian smoothing function to make the histogram continuous.  
- All descriptors are normalized by the total number of atoms in the formula.  
- Missing data in elemental properties are handled automatically.  
- The output vector can be directly used as input features for ML models (e.g., regression, classification).

---

## 5. Licensing and citation  (License, Citing)
**License(About License)**　This software is released under the MIT License, see the LICENSE.

**Citation(Citing)**  R. Jalem, M. Nakayama, Y. Noda, T. Le, I. Takeuchi, Y. Tateyama, H. Yamasaki, "A general representation scheme for crystalline solids based on Voronoi-tessellation real feature values and atomic property data", Sci. Technol. Adv. Mater., 19, 231-242 (2018) [DOI: 10.1080/14686996.2018.1439253](https://doi.org/10.1080/14686996.2018.1439253)

## Funding
Kakenhi 19H05815, 20H02436, Japan
