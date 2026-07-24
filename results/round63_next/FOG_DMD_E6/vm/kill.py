import subprocess
subprocess.run("pkill -9 -f run_e6.py; sleep 1; "
               "rm -f /content/E6_s1.json /content/E6_s2.json /content/e6_s1.log /content/e6_s2.log; "
               "echo KILLED", shell=True, executable="/bin/bash")
print("KILL_DONE")
