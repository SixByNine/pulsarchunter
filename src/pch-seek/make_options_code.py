#!/usr/bin/python


file=open("options.file")

expect="text"
short=""
long=""
argt="no_argument"
arg=""
hidden=0
code=list()
desc=list()
options=list()
text=list()
for line in file:
	if line.strip()=="" or line.strip().startswith("!!"):
		continue
	if expect=="text":
		if line.startswith("BEGIN_OPTIONS"):
			expect="option"
			continue
		text.append(line)
		continue
	
	if expect=="option":
		elems=line.split();
		if elems[0]=="*":
			elems=elems[1:]
			hidden=1
		for str in elems:
			if str.startswith("--"):
				long=str[2:]
				continue
			if str.startswith("-"):
				short=str[1]
				continue
			if str.startswith("{"):
				argt="optional_argument"
				arg=str[1:-1]
				continue
			if str.startswith("["):
				argt="required_argument"
				arg=str[1:-1]
				continue



		expect="desc"
		continue
	if expect=="desc":
		if line[0]=='{':
			code.append("\t\t"+line)
			expect="code"
			continue
		desc.append(line)
		continue
	if expect=="code":
		code.append("\t\t"+line)
		if line[0]=="}":
			options.append(dict(s=short, l=long,d=desc, c=code,a=arg,at=argt,h=hidden))

			short=""
			long=""
			arg=""
			hidden=0
			argt="no_argument"
			desc=list()
			code=list()
			expect="option"
		continue

counter=0
optargs=dict()
help=list()
help_hidden=list()
shortcode=list()
longcode=list()
shortopts=list()
for option in options:
	helpstr=""
	colons=""
	if option['h']==1:
		helpstr="*"
	if option['s']!="":
		helpstr+=" -"+option['s']
	if option['l']!="":
		helpstr+=" --"+option['l']

	if option['at']=="required_argument":
		helpstr+=("["+option['a']+"]")
		colons=":"

	if option['at']=="optional_argument":
		helpstr+=("{"+option['a']+"}")
		colons="::"



	while len(helpstr) < 20:
		helpstr+=" "
	for line in option['d']:
		helpstr+="      "
		helpstr+=line

	help_hidden.append(helpstr)
	if option['h']==0:
		help.append(helpstr)

	codestr=""
	optstr='\tlong_opt[%i].name="%s";\n'%(counter,option['l'])
	optstr+='\tlong_opt[%i].has_arg=%s;\n'%(counter,option['at'])
	if  option['s']=="":
		optstr+='\tlong_opt[%i].flag=opt_flag;\n'%(counter)
		optstr+='\tlong_opt[%i].val=0x%08x;\n'%(counter,counter)
		codestr+="\tcase 0x%08x:\n"%(counter);
		key="z%d"%counter
	else:
		optstr+='\tlong_opt[%i].flag=NULL;\n'%(counter)
		optstr+="\tlong_opt[%i].val='%s';\n"%(counter,option['s'])
		codestr+="\t\tcase '%s':\n"%option['s']
		shortopts.append(option['s']+colons)
		key=option['s']

	
	optargs[key]=(optstr)

	for line in option['c']:
		codestr+=line
	codestr+="\t\tbreak;\n"

	if  option['s']=="":
		longcode.append(codestr)
	else:
		shortcode.append(codestr)


	counter+=1


file = open("pch-seek-options.C","w")

file.write("#include <config.h>\n")
file.write("#include <stdlib.h>\n")
file.write("#include <stdio.h>\n")
file.write("#include <getopt.h>\n")
file.write("void pch_seek_help() {\n")
for str in text:
	file.write('\tprintf("%s");\n'%(str.replace("\n","\\n")))
for str in help:
	file.write('\tprintf("%s");\n'%(str.replace("\n","\\n")))
file.write("}\n\n")
file.write("void pch_seek_help_hidden() {\n")
for str in text:
	file.write('\tprintf("%s");\n'%(str.replace("\n","\\n")))

file.write('\tprintf("WARNING: Some of these options don\'t work well at all\\n");\n')
file.write('\tprintf("         Please don\'t complain if they just segfault or something!\\n");\n')
for str in help_hidden:
	file.write('\tprintf("%s");\n'%(str.replace("\n","\\n")))
file.write("}\n\n")


file.write("void pch_seek_optargs(int* nopt,struct option* long_opt, int* opt_flag){\n")
file.write("\t*nopt=%d;\n"%len(optargs))
keys=optargs.keys()
keys.sort()
for key in keys:
	str=optargs[key]
	file.write(str)
file.write("}\n")

file.close()


file = open("pch-seek-options.h","w")

shortopts.sort()
file.write('const char* args="')
for o in shortopts:
	file.write(o)
file.write('";\n')
file.write("char c;\nint nopt, opt_flag;\n")
file.write("struct option long_opt[%d];"%len(optargs))
file.write("pch_seek_optargs(&nopt,long_opt,&opt_flag);\n")
file.write("while ((c = getopt_long(argc, argv, args, long_opt, &nopt)) != -1) {\n")
file.write("\t// handle args\n")
file.write("\tswitch (opt_flag) {\n")
for str in longcode:
	 file.write(str)
file.write("\tdefault:\n")
file.write("\t\tswitch(c){\n")
for str in shortcode:
	file.write(str)
file.write("\t\t}\n")
file.write("\t};\n")
file.write("opt_flag=0;\n")
file.write("}\n")
file.close()

	
