# Degradation-aware isoform quantification

daiso is tool for isoform abundance estimation from ONT direct RNA-seq that models degradation bias. 

## Installation

Clone the package:
```
git clone https://github.com/jleechung/daiso
cd daiso
```

Set up a python virtual environment:
```
python3 -m venv myenv
source ./myenv/bin/activate
```

Install requirements:
```
pip3 install -r requirements.txt
```

Test the installation:
```
python3 daiso.py -h
```
This should print the help message.

## Usage

