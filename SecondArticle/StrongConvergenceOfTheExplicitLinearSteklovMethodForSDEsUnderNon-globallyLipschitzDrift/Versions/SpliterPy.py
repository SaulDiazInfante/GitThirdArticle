f = open('StrongConvergenceLSMethod.tex', 'r')
strSectionName = 'foo.tex'
sec = open(strSectionName, 'w')
count = 0
gi = open('SecondArticle.tex', 'w')
swsc = False
swgi = True
for line in f:
	count = count + 1
	if line.find('\section{') == 0:
		sec.close()
		swgi = False
		strSectionName = line[line.index('{')+1:line.index('}')]
		strSectionName = strSectionName.replace(' ','')
		strSectionName = strSectionName+'.tex'
		latexin ='\t\t\\input{./'+strSectionName+'}\n'
		sec = open(strSectionName, 'w')
		gi.write('\t'+line)
		gi.write(latexin)
		
		print line
		swsc = True
		swgi = False
	
	if swgi:
		gi.write(line)
	
	if swsc:
		
		if line.find('\end{document}') == 0:
			gi.write(line)
		else:
			if line.find('\section{')==-1:
				sec.write(line)

gi.close()
sec.close()
f.close()
print 'Line Count:', count

