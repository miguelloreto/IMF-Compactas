import numpy as np
import pandas as pd
import os

out=os.listdir(os.getcwd())

for i in out:
    i=i.split("-")
    if(i[0] == "spec"):
        text=open('-'.join(i), 'r')
        linhas=text.readlines()
        head = linhas[39]
        head = head.split("j")
        if(head[0]!= '   # '):
            print('-'.join(i))
