from setuptools import setup
import os

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='bioseq',  
     version='1.0',
     author="antonioalmeida",
     description="Handle biological sequences and perform common Bioinformatics operations",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/antonioalmeida/fcup-abi",
     packages=setuptools.find_packages(),
     install_requires=required,
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )