import subprocess
cmd = ("cd /content && setsid nohup python run_e6.py --arms A4,A0,A1 "
       "--seeds 0,1,2,3,4 --out /content/E6_s1.json --dev cuda "
       "> /content/e6_s1.log 2>&1 & echo $! > /content/e6_s1.pid")
subprocess.run(cmd, shell=True, executable="/bin/bash")
import time; time.sleep(2)
print("LAUNCHED_S1 pid", open("/content/e6_s1.pid").read().strip())
