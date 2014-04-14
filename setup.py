#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
import os

install_requires = ['numpy>=1.3.0',]

setup(
    name = "dmgraphanalysis",
    version = '0.0.1dev',
    packages = ['dmgraphanalysis'],
    install_requires=install_requires,
    author = "David Meunier",
    description = "Graph analysis for nipype"
)

