# CSA S16 Python

**CSA_S16_python** is a Python package for structural engineers working with CSA S16:24. It includes a comprehensive set of functions to calculate factored resistances member and connection checks — all in accordance with CSA S16.

This package is designed for use in Jupyter notebooks and web applications, with support for rendered equations via the `handcalcs` and `IPython.display` libraries.

---

## Features

- ...

---

## Installation

```bash
pip install CSA-S16-python==0.1.0
```
---

## Use

```python
from CSA_S16 import *
help(CSA_S16) # see available functions
```

Example of Use
```python
T_r = T_r_y_func(1000 * mm, 350 * MPa)
display(Math(T_r[0])) # see LaTeX from first item
```
```python
T_r[1] # get numerica value from second item
```
