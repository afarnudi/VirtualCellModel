#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    subprocess.call('cd ..; doxygen', shell=True)

import sphinx_rtd_theme

html_theme = "sphinx_rtd_theme"

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

def setup(app):
    app.add_stylesheet("main_stylesheet.css")

extensions = ['breathe']
breathe_projects = { 'The Virtual Cell Model': '../xml' }
breathe_default_project = "The Virtual Cell Model"
templates_path = ['_templates']
html_static_path = ['_static']
source_suffix = '.rst'
master_doc = 'index'
project = 'The Virtual Cell Model'
copyright = '2016, Mohammad Reza Ejtehadi, Morteza Mahmoudi, Maziar Heidari,Tiam Heydari, Ali Farnudi'
author = 'Ali Farnudi'

html_logo = 'logo_1.png'

exclude_patterns = []
highlight_language = 'c++'
pygments_style = 'sphinx'
todo_include_todos = True
htmlhelp_basename = 'TVCMdoc'
