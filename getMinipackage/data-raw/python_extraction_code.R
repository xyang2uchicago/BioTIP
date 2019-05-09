##
## #!/bin/python
##
##
## # Replace PATH_FILE your local file directory
## gtf = ('PATH_FILE.gtf')
## outF = open('GENCODEv25_gtf_summary_21.txt','w')
##
## def getquote(str,f,target):
##     targetLen = len(target)+2
##     strInd = str.find(target)
##     st = strInd + len(target)+2
##     ed = st + str[st:].find('";')
##     #print(st,ed)
##     f.write('\t'+str[st:ed]) if strInd!= -1 else f.write('\t'+'NA.')
##
## with open(gtf, 'r') as f:
##    for line in f:
##    if line[0] != '#':
##    chromosome = line.split('\t')[0]
##    st = line.split('\t')[3]
##    ed = line.split('\t')[4]
##    strand = line.split('\t')[6]
##    type = line.split('\t')[2]
##    outF.write(chromosome+'\t'+st+'\t'+ed+'\t'+strand+'\t'+type)
##    c = 'transcript_id'
##    g = 'gene_name'
##    t = 'transcript_type'
##    getquote(line,outF,c)
##    getquote(line,outF,g)
##    getquote(line,outF,t)
##    outF.write('\n')
## outF.close()

