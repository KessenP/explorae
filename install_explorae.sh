#!/bin/bash

set -e

# Installation des dépendances Python
echo "Installation des dépendances du requirements.txt..."
pip install --break-system-packages -r requirements.txt

echo "Installation de pyrosetta-installer..."
pip install --break-system-packages pyrosetta-installer

echo "Installation de PyRosetta..."
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

echo "Installation terminée !"