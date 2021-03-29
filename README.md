# BOSS SBI 
[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/changhoonhahn/molino/blob/main/LICENSE)
[![Gitter](https://badges.gitter.im/boss_sbi/community.svg)](https://gitter.im/boss_sbi/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

Simulation-based Inference of the BOSS survey


Details and resources on the original BOSS analysis: [https://sites.google.com/view/learningouruniverse/boss-analysis?authuser=0](https://sites.google.com/view/learningouruniverse/boss-analysis?authuser=0)



## Installation 
First [set up a new anaconda environment](#setting-up-a-conda-environment) to
avoid any package conflicts. 

Activate the conda environment and then install the `boss_sbi` package
```bash
# first clone the repo
git clone https://github.com/changhoonhahn/boss_sbi.git

# go to the repo
cd boss_sbi

# install the package 
pip install -e . 
```

You will also need to set the `$QUIJOTE_DIR` environment variable. To do this, add the following line 
```
export QUIJOTE_DIR="/projects/QUIJOTE/"
```
to your `~/.bashrc` file. If you don't know how to do this then just copy paste the following: 
```
echo 'export QUIJOTE_DIR="/projects/QUIJOTE/"' >> ~/.bashrc
```

Once you've added the line, don't forget to run
```
source ~/.bashrc
```

### Setting Up a Conda Environment 
#### On `tiger` 
If you're on Princeton's `tiger` cluster, you don't have to install anaconda.
You can load it using 
```bash
module load anaconda 
```

Afterwards you can create a new conda environment using
```bash
conda create -n ENV_NAME_HERE python=3.7 ipython 
```
and following the instructions. 


To activate the conda environment you created
```
conda activate ENV_NAME_HERE 
```

Later, if you want to exist the conda environemtn
```bash
conda deactivate 
```


### Dependencies
The `boss_sbi` package requires the following python pacakges: 
- [nbodykit](https://nbodykit.readthedocs.io/) 
- [pymangle](https://github.com/esheldon/pymangle)

**tl;dr** Run the following lines after activating the conda environment 
```
conda install -c bccp nbodykit
pip install pymangle
```


## Generating an HOD catalog for HR Quijote LHC 
```python
import numpy as np 
from boss_sbi.halos import Quijote_LHC_HR
from boss_sbi import galaxies as Galaxies

# read in halo catalog 
halos = Quijote_LHC_HR(1, z=0.5)

# get LOWZ HOD parameters
theta_hod = Galaxies.thetahod_lowz_ngc()

# apply HOD 
gals = Galaxies.hodGalaxies(halos, theta_hod, seed=0) 
print(np.array(gals['Position']))
```
