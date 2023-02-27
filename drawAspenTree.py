from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, AttrFace, TextFace, faces
from optparse import OptionParser
from palettable.colorbrewer.qualitative import Paired_12

## example xvfb-run python drawTree.py -t data/little_ksu/h3.l.r.newick -o output/treetestlk.svg -e 
## output/little_ksu/nsyn/skip_stoppers/trees/etetest/h3_treescheme_291_Node_3506
parser = OptionParser()
parser.add_option("-t", "--treefile", dest="treefile",
                  help="newick tree file")
parser.add_option("-e", "--eventfile", dest="eventfile",
                  help="events file")
parser.add_option("-p", "--passagefile", dest="passagefile",help="tab delimited file without a header: Nodename\t[0123]")
parser.add_option("-o", "--output", dest="output",
                  help="output file")
parser.add_option("-s", "--scale", dest="scale", type="int",
                  help="8000 for little ksu, 3 for krya")
parser.add_option("-c", "--circle_size", dest="size", type="float",
                  help="20 for little ksu, ? for krya")
parser.add_option("-w", "--width", dest="width", type="int",
                  help="870 for little ksu, ? for krya")
parser.add_option("-l", "--large", dest="large",action="store_true", default=False,
                  help="do not set vertical margin between branches( branch margins used for scheme)")
parser.add_option("-n", "--notext", dest="notext",action="store_true", default=False, help="do not add derived codon labels to internal nodes")
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help="print node names")
(options, args) = parser.parse_args()
t = Tree(options.treefile, format=1) #flexible with internal node names
farthest, dist = t.get_farthest_node()
print ("The farthest node from root is", farthest.name, "with dist=", dist)
if not options.scale:
    options.scale = int(1700/dist) #900/dist
#root = t.get_tree_root()
#farthest = root.get_farthest_leaf()
#print "The height is ", t.get_distance(root.name, farthest.name)

# Basic tree style
ts = TreeStyle()
ts.show_leaf_name = False
if options.large:
   ts.show_scale = False
#ts.optimal_scale_level = "full"
ts.scale = options.scale #3 for krya # pixels per branch length unit
if not options.large:
    ts.branch_vertical_margin = 2


passagecolor = {'0':"blue",'1':"yellow",'2':"red",'3':"grey"}
nodepassage = {}
if options.passagefile:
	print("Passagefile "+options.passagefile)
	pf = open(options.passagefile, "r")
	pfc = pf.readlines()
	for line in pfc:
		nname = line.split('\t')[0]
		state = line.split('\t')[1].rstrip('\n')
		nodepassage[nname] = state
print(nodepassage)
efile = open(options.eventfile, "r")
subtrees = []
syns = {}
derived = {}
counts = {}
rootstring = efile.readline()
rootnode = rootstring.split(':')[1].split('|')[0]
derived[rootnode] = rootstring.split('|')[1].rstrip('\n')
roottext = rootstring.split('|')[2].rstrip('\n')
e = {"anc":rootstring.split(':')[1].rstrip('\n'),"events":[],"subtree":[],"text":roottext}
subtrees.append(e)
while 1:
	line = efile.readline()
	if not line:
		break
	if line[0:3] == "Syn":
		synsubs = line.split(':')[1].rstrip('\n,').split(',')
		for ss in synsubs:
			if ss:
				syns[ss.split('|')[0]] = ss.split('|')[1]
		break
	anc = line.split(':')[1].rstrip('\n')
	ancletter = anc.split('|')[1]
	anctext = anc.split('|')[2]
	tevents = efile.readline().rstrip('\n,').split(':')[1].split(',')
	events = []
	for e in tevents:
		#print(e)
		nod = e.split('|')[0]
		events.append(nod)
		letter = e.split('|')[1]
		derived[nod] = letter
		if letter in counts:
			counts[letter] = counts[letter]+1
		else:
			counts[letter] = 1
		
	subtree = efile.readline().rstrip('\n,').split(':')[1].split(',')
	#print(ancletter)
	if ancletter in counts:
		counts[ancletter] = counts[ancletter]+len(subtree)
	else:
		counts[ancletter] = len(subtree)
	e = {"anc":anc,"events":events,"subtree":subtree,"text":anctext}	
	subtrees.append(e)
	
dercolor = {}
i = 0

for der, c in sorted(counts.items(),  key=lambda kv: (-kv[1], kv[0])):
		dercolor[der] = Paired_12.hex_colors[i%len(Paired_12.hex_colors)]
		#print(der)
		#print(c)
		#print(dercolor[der])
		i = i+1
dercolor[derived[rootnode]] = "lightgray"
#print("root aa" + derived[rootnode])
#for der in set(derived.values()):
#		dercolor[der] = Paired_12.hex_colors[i%len(Paired_12.hex_colors)]
#		print(der)
#		print(dercolor[der])
#		i = i+1

def splitNode(n):
        parent = n.up
        length = n.dist       
        detached = n.detach()
        newnode = parent.add_child(name = detached.name)
        newnode.dist = length/2
        detached.name = detached.name+"child"
        grafted = newnode.add_child(detached)
        grafted.dist = length/2
        return grafted 

		
## todo : traverse and color grey
		
for e in subtrees:		
	if not options.large:
		eventnodes = []
		allnodes = []
		for n in t.traverse():
		  if not n.is_root():
				  allnodes.append(n)  
		  if n.name in e["events"]:
				  eventnodes.append(n)
		for n in allnodes:
		  nname = n.name
		  newchild = splitNode(n)
		  if nname in e["subtree"] and not nname in e["events"]:
				  e["subtree"].append(newchild.name) 




# Creates an independent node style for each node, which is
# initialized with a red foreground color.
t.ladderize()

prunedSyns = {}
for snode in syns.keys():
		found = 0
		for e in subtrees:
			nname = e["anc"].split('|')[0]
			n = t.search_nodes(name=nname)[0]
			while n:
				if n.name == snode:
					found = 1
					prunedSyns[snode] = syns[snode]
					break
				n = n.up
			if found == 1:
				break			
#print(syns)
#print(prunedSyns)
#treecolor = "lightseagreen" #"mediumseagreen" # "lightseagreen" # "mediumaquamarine" #"seagreen" "lightseagreen"
backcolor = dercolor[derived[rootnode]] #"lightgray" # "silver" #"DarkGray"
#eventscolor = "#006060" #"#707070" # "#006060"  #"teal" "Black"
for n in t.traverse():
	nstyle = NodeStyle()
	nstyle["hz_line_color"] = backcolor
	nstyle["vt_line_color"] = backcolor
	nstyle["size"] = 0
	nstyle["vt_line_width"] = 1
	nstyle["hz_line_width"] = 1
	for e in subtrees:
		treecolor = dercolor[e["anc"].split('|')[1]]
		#print(treecolor)
		if n.name != e["anc"].split('|')[0]:
			#print ("not an ancestor")
			if n.name in e["events"]:
				eventscolor = dercolor[derived[n.name]]	
				nstyle["fgcolor"] = eventscolor 
				nstyle["size"] = options.size
				nstyle["hz_line_color"] = treecolor #
				nstyle["vt_line_color"] = treecolor #backcolor
				print("parsing events..")
				if options.passagefile:
					print ("Will add passage labels")
				#	n.img_style["bgcolor"] = passagecolor[nodepassage[n.name]]
				#	n.img_style["size"] = 12
					nstyle["bgcolor"] = passagecolor[nodepassage[n.name]]
				if (options.debug):
					print ("Will add node labels")
					n.add_face(TextFace(n.name), column=0)
			if n.name in e["subtree"]:
				nstyle["hz_line_color"] = treecolor
				nstyle["vt_line_color"] = treecolor
		else:
			eventscolor = dercolor[derived[n.name]]	
			nstyle["fgcolor"] = eventscolor
			nstyle["size"] = options.size
			if (options.debug):
					print ("Will add node labels")
					n.add_face(TextFace(n.name), column=0)
			if options.passagefile:
					print ("Will add passage labels")
					nstyle["bgcolor"] = passagecolor[nodepassage[n.name]]							
			if (len(e["events"])>1 or n.name == rootnode) and len(e["anc"].split('|')[1]) == 1 and not options.notext:
				n.add_face(TextFace(e["text"]), column=1)
	n.set_style(nstyle)
	if n.name in prunedSyns.keys():
		n.add_face(TextFace(prunedSyns[n.name]), column=1)

for aa in dercolor.keys():
	ts.legend.add_face(CircleFace(10, dercolor[aa]), column=0)
	ts.legend.add_face(TextFace(aa), column=1)
	
t.render(options.output + ".pdf", tree_style=ts, dpi=300, w=options.width, units="mm")
t.render(options.output + ".png", tree_style=ts, dpi=300, w=options.width, units="mm")
t.render(options.output + ".svg", tree_style=ts, dpi=300, w=options.width, units="mm")
