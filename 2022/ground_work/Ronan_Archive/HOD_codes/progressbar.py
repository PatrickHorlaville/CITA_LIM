import sys
import numpy as np

def progress_bar(i, N, bar_length):

    itick = N / bar_length
    tick = float(N) / bar_length
    s_len = int(np.round((i+1) / tick))
    s_per = int(np.round(100. * float(i+1) / N))
    s_pro = ( " |" + (s_len*"=").ljust(bar_length) + "|" + 
              str(s_per) + "%" + 
              " (" + str(i+1) + "/" + str(N) + ")" + "\r" )
    if (i%itick==0 or (i+1)==N or i==0):
        sys.stdout.flush()
        sys.stdout.write(s_pro)

    if (i+1) == N: print("")



