# VM-side launcher for cbt_camp2: JOB 3 CBT_PRIVACY_TEST (conditional-MI certification).
# Places files, verifies md5s, background-runs the privacy certification (B=16384, N=6000).
import hashlib
import os
import subprocess

BASE = '/content/results/round63_next'
DIRS = {
    'CURVED_BLINDNESS_TEST': 'curved_blindness_test.py',
    'SCRAMBLE_EXT': 'scramble_toy.py',
    'JET_TEST': 'jet_test.py',
    'CBT_LAW_PLANE': 'cbt_law_plane.py',
    'CBT_PRIVACY_TEST': 'cbt_privacy_test.py',
}
for d, fn in DIRS.items():
    os.makedirs(os.path.join(BASE, d), exist_ok=True)
    src = os.path.join('/content', fn)
    if os.path.exists(src):
        dst = os.path.join(BASE, d, fn)
        with open(src, 'rb') as a, open(dst, 'wb') as b:
            b.write(a.read())


def md5(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


for d, fn in DIRS.items():
    p = os.path.join(BASE, d, fn)
    print('MD5 %-24s %s' % (fn, md5(p) if os.path.exists(p) else 'MISSING'))
try:
    import torch
    print('TORCH', torch.__version__, 'CUDA', torch.cuda.is_available(),
          torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')
except Exception as e:
    print('env probe error:', e)

wd = os.path.join(BASE, 'CBT_PRIVACY_TEST')
log = os.path.join(wd, 'camp2_progress.log')
sh = '/content/run_camp2.sh'
with open(sh, 'w') as f:
    f.write('#!/bin/bash\ncd %s\n' % wd)
    f.write('echo "=== camp2 start $(date -u +%%H:%%M:%%S) ===" >> %s\n' % log)
    f.write('python -u cbt_privacy_test.py --mode run >> %s 2>&1\n' % log)
    f.write('echo "CAMP2_ALL_DONE $(date -u +%%H:%%M:%%S)" >> %s\n' % log)
open(log, 'w').close()
p = subprocess.Popen(['bash', sh], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
print('LAUNCHED_CAMP2_PID', p.pid)
