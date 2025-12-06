# twoface: transit lightcurves for asymmetric planets transiting spotted stars

This repository contains `twoface`, an exoplanet transit model for asymmetric planets transiting spotted stars. The code is named after the Batman villain Two-Face, reflecting its ability to model planets with asymmetric, two-sided silhouettes as they cross spotted, troublesome stellar surfaces. It builds upon the numerical integration framework and architecture of `spotrod` (`Béky, Kipping & Holman 2014 <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.3686B/abstract>`_) and the geometry of `catwoman` (`Jones & Espinoza 2020 <https://doi.org/10.21105/joss.02382>`_, `Espinoza & Jones 2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210615687E/abstract>`_).


## Installation

Install the latest version of `twoface` using `pip`:

```
pip install git+https://github.com/mans-holmberg/twoface
```

## Usage

```python
import matplotlib.pyplot as plt
import numpy as np
import twoface

t = np.linspace(-0.8, 0.8, 1000)
params = twoface.TransitParams()

params.per = 3.
params.t0 = 0.
params.a = 1.5
params.inc = 90.
params.ecc = 0.
params.w = 90.
params.rp = 0.1
params.rp2 = 0.08
params.phi = -45.

params.limb_dark = 'quadratic'
params.u = [0.3, 0.2]

params.num_spots = 2
params.spot_mu = [0.99, 0.7]
params.spot_phi = [1.6, 0.1]
params.spot_radius = [0.1, 0.1]
params.spot_f = [0.8, 0.8]

model = twoface.TransitModel(params, t)
f = model.light_curve(params)

plt.plot(t, f)
plt.xlabel("Time from central transit/days")
plt.ylabel("Relative flux")
plt.show()
```