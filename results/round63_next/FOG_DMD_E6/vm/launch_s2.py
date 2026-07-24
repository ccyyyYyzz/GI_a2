import subprocess
cmd = ("cd /content && setsid nohup python run_e6.py --arms A2,A3 "
       "--seeds 0,1,2,3,4 --out /content/E6_s2.json --dev cuda "
       "> /content/e6_s2.log 2>&1 & echo $! > /content/e6_s2.pid")
subprocess.run(cmd, shell=True, executable="/bin/bash")
import time; time.sleep(2)
print("LAUNCHED_S2 pid", open("/content/e6_s2.pid").read().strip())
