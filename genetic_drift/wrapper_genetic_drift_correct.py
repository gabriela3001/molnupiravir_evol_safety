import subprocess
import pickle

script_path = 'script_genetic_drift_correct.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null -o log.txt python ' + script_path

with open('missing_files.txt', 'rb') as f:
  missing = pickle.load(f)

for param in range(9261):
#for param in missing:
	paramstr = ' '.join(map(str,[param]))
	_cmd = ' '.join([cmd, paramstr])
	subprocess.run(_cmd, shell=True)