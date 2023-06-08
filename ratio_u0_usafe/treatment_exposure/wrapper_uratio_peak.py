import subprocess
import numpy as np

script_path = 'script_uratio_treatment_exposure.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path


for mind in range(800):
    paramstr = ' '.join(map(str,[mind]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)