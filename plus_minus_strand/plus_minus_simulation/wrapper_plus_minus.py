import subprocess

script_path = 'script_plus_minus.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for param in range(6):
    for niter in range(1000):
        paramstr = ' '.join(map(str,[param, niter]))
        _cmd = ' '.join([cmd, paramstr])
        subprocess.run(_cmd, shell=True)