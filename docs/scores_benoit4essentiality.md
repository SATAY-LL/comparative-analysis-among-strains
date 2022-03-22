---
marp: true
theme: default
author: Leila Inigo de la Cruz
paginate: true
style: |
    img[alt~="center"] {
      display: block;
      margin: 0 auto;
    }

---

# Domain Likelihood scores 

---

## Using the idea of Benoit in: https://elifesciences.org/articles/23570

![w:700 center](scores_Benoit.png)

---

# Python implementation steps

- Import data 

![w:900 center](main-data.png)

---

- Take out the genes that have multiple genetic locations : e.g 'ECM1', 'MYO4', 'SSA1', 'RPL30', 'RPL19A', 'PRE8', 'PRE7','RPL19B', 'ELP6', 'ETR1', etc.

- Compute the scores 

![w:800 center](scores-wrokflow_1.png)

---

![w:650 center](scores-wrokflow_2.png)


---

## Implementation of the critical step towards the score


- The largest interval from tr=n to tr=n+t. Example

```python
t=5
vector_insertions=[332164, 332194, 332241, 332294, 332649, 332652, 332700]

vector_insertions[t:]= [332652, 332700]
vector_insertions[:-t]=[332164, 332194]

interval=vector_insertions[t:]-vector_insertions[:-t]=[[488, 506]]

max_interval=max(interval)=506bp

```

---

![w:650 center](fig_wt_new_pairplot_dls_with_essentials.png)