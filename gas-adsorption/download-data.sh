#!/bin/sh

set -ev

python3 isothermal-ads-gas-save-csv.py 1 10000 298 methane

python3 isothermal-ads-gas-save-csv.py 1 10000 298 hydrogen
python3 isothermal-ads-gas-save-csv.py 1 10000 77 hydrogen
python3 isothermal-ads-gas-save-csv.py 1 10000 160 hydrogen

