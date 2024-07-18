import pandas as pd
import numpy as np
import sys
import pysam
input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

for i, row in input_bed.iterrows():
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    col12 = []
    col13 = []
    col14 = []
    col15 = []
    col16 = []
    col17 = []
    col18 = []
    col19 = []
    col20 = []
    col21 = []
    col22 = []
    col23 = []
    col24 = []
    col25 = []
    col26 = []
    col27 = []
    col28 = []
    if str(row[0]) in contigs:
        for gc in input_anno.fetch(str(row[0]), row[1], row[2], parser = pysam.asTuple()):
            if gc[3] != "":
                col4.append(int(gc[3]))
            if gc[4] != "":
                col5.append(int(gc[4]))
            if gc[5] != "":
                col6.append(int(gc[5]))
            if gc[6] != "":
                col7.append(int(gc[6]))
            if gc[7] != "":
                col8.append(int(gc[7]))
            if gc[8] != "":
                col9.append(int(gc[8]))
            if gc[9] != "":
                col10.append(int(gc[9]))
            if gc[10] != "":  
                col11.append(int(gc[10]))
            if gc[11] != "":
                col12.append(int(gc[11]))
            if gc[12] != "":    
                col13.append(int(gc[12]))
            if gc[13] != "":    
                col14.append(int(gc[13]))
            if gc[14] != "":    
                col15.append(int(gc[14]))
            if gc[15] != "":    
                col16.append(int(gc[15]))
            if gc[16] != "":    
                col17.append(int(gc[16]))
            if gc[17] != "":    
                col18.append(int(gc[17]))
            if gc[18] != "":    
                col19.append(int(gc[18]))
            if gc[19] != "":    
                col20.append(int(gc[19]))
            if gc[20] != "":    
                col21.append(int(gc[20]))
            if gc[21] != "":    
                col22.append(int(gc[21]))
            if gc[22] != "":    
                col23.append(int(gc[22]))
            if gc[23] != "":    
                col24.append(int(gc[23]))
            if gc[24] != "":   
                col25.append(int(gc[24]))
            if gc[25] != "":    
                col26.append(int(gc[25]))
            if gc[26] != "":    
                col27.append(int(gc[26]))
            if gc[27] != "":    
                col28.append(int(gc[27]))
        if col4 == []:
            print(row[0], row[1], row[2], ".", ".", ".", ".", ".", ".",".", ".", ".", ".", ".", ".", ".", ".", ".", ".",".", ".", ".", ".", ".", ".", ".", ".", ".", sep="\t")
        else:
            print(row[0], row[1], row[2],
                    max(col4), max(col5), max(col6),
                    max(col7), max(col8), max(col9),
                    max(col10), max(col11), max(col12), 
                    max(col13), max(col14), max(col15), max(col16),
                    max(col17), max(col18), max(col19), max(col20),
                    max(col21), max(col22), max(col23), max(col24),
                    max(col25), max(col26), max(col27), max(col28), sep="\t")
    else:
        print(row[0], row[1], row[2], ".", ".", ".", ".", ".", ".",".", ".", ".", ".", ".", ".", ".", ".", ".", ".",".", ".", ".", ".", ".", ".", ".", ".", ".", sep="\t")