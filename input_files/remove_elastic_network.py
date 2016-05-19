#!/usr/bin/python


#######################################################################
## This script is to generate itp file remove elastic network #########
### Usage: python remove_elastic_network.py input_filename.itp output_filename.itp ###
#######################################################################

import sys

input_ITP_file_name = sys.argv[1]
input_ITP_file = open(input_ITP_file_name, 'r')

#output_ITP_file_name = sys.argv[2]
output_ITP_file = open(sys.argv[2], 'w')

output_removed_lines_file=open("output_removed_lines.dat", "w")
aa_name = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HIE","HIP","HID","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]

#sec1 = [1:1085]
#sec2 = [1086:2452]


for line1 in input_ITP_file.readlines():
    
    line=line1.split()
    aa= line1[37:40].strip()
    ba= line1[40:43].strip()
    da= line1[52:55].strip()

    print(line)
    print(line1)
    #sec1_an = int(line.split()[0])
    #sec2_an = int(line.split()[1])
    
#    68-92

    if (line1.find("RUBBER_FC*1.000000")>0) and ((int(line[0]) > 78 and int(line[0]) < 95) or (int(line[1]) > 78 and int(line[1]) < 95)): 
        output_removed_lines_file.write(line1)

    else:
        output_ITP_file.write(line1)
       # output_removed_lines[i] = remove_line
       # if line.startswith("#endif"):
  
input_ITP_file.close()
output_ITP_file.close()
output_removed_lines_file.close()


#open(output_removed_lines_file, 'w').writelines(output_removed_lines)
