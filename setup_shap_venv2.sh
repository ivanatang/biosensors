#!/bin/bash
set -e
module load anaconda
VENV_DIR="/projects/ivta1597/software/shap_venv"
rm -rf "$VENV_DIR"
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
pip install --upgrade pip
pip install shap scikit-learn pandas numpy matplotlib scipy openpyxl
echo "=== shap venv versions ==="
python3 -c "
import numpy, pandas, sklearn, shap, matplotlib, scipy, openpyxl
print('numpy:', numpy.__version__)
print('pandas:', pandas.__version__)
print('sklearn:', sklearn.__version__)
print('shap:', shap.__version__)
print('matplotlib:', matplotlib.__version__)
print('scipy:', scipy.__version__)
print('openpyxl:', openpyxl.__version__)
"
echo "DONE_MARKER_OK"
