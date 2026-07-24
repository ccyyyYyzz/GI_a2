# VM-side launcher: place the 3 byte-identical files into the frozen directory structure,
# verify md5s, then launch curved_blindness_test.py in the BACKGROUND on the A100 so the CLI
# exec returns immediately and we can poll run.log.  NO edits to the frozen script.
import os, subprocess, shutil, hashlib
BASE = '/content/results/round63_next'
cbt = os.path.join(BASE, 'CURVED_BLINDNESS_TEST')
scr = os.path.join(BASE, 'SCRAMBLE_EXT')
jet = os.path.join(BASE, 'JET_TEST')
for d in (cbt, scr, jet):
    os.makedirs(d, exist_ok=True)
shutil.copy('/content/curved_blindness_test.py', os.path.join(cbt, 'curved_blindness_test.py'))
shutil.copy('/content/scramble_toy.py', os.path.join(scr, 'scramble_toy.py'))
shutil.copy('/content/jet_test.py', os.path.join(jet, 'jet_test.py'))


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


print('MD5 cbt      =', md5(os.path.join(cbt, 'curved_blindness_test.py')))
print('MD5 scramble =', md5(os.path.join(scr, 'scramble_toy.py')))
print('MD5 jet      =', md5(os.path.join(jet, 'jet_test.py')))
try:
    import torch
    print('TORCH', torch.__version__, 'CUDA', torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')
    import numpy, scipy
    print('numpy', numpy.__version__, 'scipy', scipy.__version__)
except Exception as e:
    print('env probe error:', e)
logf = open(os.path.join(cbt, 'colab_run.out'), 'w')
p = subprocess.Popen(['python', '-u', 'curved_blindness_test.py'], cwd=cbt,
                     stdout=logf, stderr=subprocess.STDOUT)
print('LAUNCHED_PID', p.pid)
