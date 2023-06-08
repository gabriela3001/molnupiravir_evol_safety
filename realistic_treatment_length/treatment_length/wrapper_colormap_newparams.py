import subprocess

script_path = 'script_colormap_newparams.py'
cmd = 'bsub -q new-short -R "rusage[mem=1000]" -u /dev/null python ' + script_path

for param in [6836,6837,6838,6839]:
    paramstr = ' '.join(map(str,[param]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)