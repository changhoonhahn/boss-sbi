# BOSS SBI 
[![MIT License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/changhoonhahn/molino/blob/main/LICENSE)

Simulation-based Inference of the BOSS survey


Details and resources on the original BOSS analysis: [https://sites.google.com/view/learningouruniverse/boss-analysis?authuser=0](https://sites.google.com/view/learningouruniverse/boss-analysis?authuser=0)



## Installation 
To install the `boss_sbi` package
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
