import subprocess

script_path = 'script_lethal_defection.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for param in range(9261):
	paramstr = ' '.join(map(str,[param]))
	_cmd = ' '.join([cmd, paramstr])
	subprocess.run(_cmd, shell=True)