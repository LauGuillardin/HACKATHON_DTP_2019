#!/usr/bin/python

i = open("SRR5218242.1.fastq")
o = open("output_1.fa", "w")
p = open("output_2.fa", "w")

count = 0 

for line in i:
	if count%4 ==0:
		if "/1" in line:
			o.write(line)
			set==1
		if "/2" in line:
			p.write(line)
			set==2
	else:
		if set ==1:
			o.write(line)
		if set ==2:
			p.write(line)
	count = count + 1
