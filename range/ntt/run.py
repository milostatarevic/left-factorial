import subprocess, sys, os, math
from time import sleep


max_cores = 12


def getmem(pid):
    """ returns current memory used by given pid in megabytes """
    proc = subprocess.Popen("ps -p %d -o rss" % pid, shell=True, stdout=subprocess.PIPE)
    return float(proc.communicate()[0].split()[-1]) / 1024.0



def run(target, size, cores):
    """ runs target == 'gmp' or 'ntt' with given # cores for size megabytes
    returns (time in seconds, maximum memory in megabytes) """
    n = int(1048576.0 * size / 8)    # size in limbs
    if target == "gmp":
        command = "./test-gmp gmp"
    elif target == "mpir":
        command = "./test-mpir gmp"
    elif target == "ntt":
        command = "./test-gmp ntt"
    else:
        assert(False)
        
    command = "%s %d %d" % (command, n, cores)
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    pid = proc.pid
    mem = 0
    while proc.poll() is None:
        mem = max(mem, getmem(proc.pid))
        sleep(0.1)
    time = float(proc.communicate()[0].split()[-1])
    return (time, mem)
    


ratio = (1 + math.sqrt(5)) / 2

size = 10.0     # input size in megabytes

format = "  %4s       %2d     %8.1f      %6.1f       %5d"
print "target    # cores    size (MB)    time (s)    mem (MB)"
print

while size < 1024.0:
    n = int(1048576 * size / 8)    # size in limbs

    L = [max_cores]
    while L[-1] > 1:
        L.append(int(math.ceil(L[-1] / 2.0)))
    regimes = [("gmp", 1), ("mpir", 1)] + [("ntt", cores) for cores in reversed(L)]

    for (target, cores) in regimes:
        time, mem = run(target, size, cores)
        print format % (target, cores, size, time, mem)

    size *= ratio
    print
    
