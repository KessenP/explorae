#!/bin/bash

set -e

# Installation des dépendances Python
echo "Installation des dépendances du requirements.txt..."
pip install -r requirements.txt

echo "Installation de pyrosetta-installer..."
pip install pyrosetta-installer

echo "Installation de PyRosetta..."
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

echo "END"