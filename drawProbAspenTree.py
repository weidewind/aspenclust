from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace
from optparse import OptionParser
from palettable.colorbrewer.qualitative import Paired_12
import re

## example xvfb-run python drawTree.py -t data/little_ksu/h3.l.r.newick -o output/treetestlk.svg -e 
## output/little_ksu/nsyn/skip_stoppers/trees/etetest/h3_treescheme_291_Node_3506
parser = OptionParser()
parser.add_option("-t", "--treefile", dest="treefile",
                  help="newick tree file")
parser.add_option("-e", "--eventfile", dest="eventfile",
                  help="events file")
parser.add_option("-o", "--output", dest="output",
                  help="output file")
parser.add_option("-p", "--prune", dest="prune", type="float",default=0,help="Do not draw subs with probability less than this value")
parser.add_option("-s", "--scale", dest="scale", type="int",
                  help="8000 for little ksu, 3 for krya")
parser.add_option("-c", "--circle_size", dest="size", type="float",
                  help="20 for little ksu, ? for krya")
parser.add_option("-w", "--width", dest="width", type="int",
                  help="870 for little ksu, ? for krya")
parser.add_option("-l", "--large", dest="large",action="store_true", default=False,
                  help="do not set vertical margin between branches( branch margins used for scheme)")
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

efile = open(options.eventfile, "r")
events = {}
counts = {}

# line: Node_2344:ATT(F)->M|0.132,ATC(F)->M|0.342,
while 1:
	line = efile.readline()
	if not line:
		break
	nname = line.split(':')[0].rstrip('\n')
	tevents = line.split(':')[1].rstrip('\n,').split(',')
	for e in tevents:
		print(e)
		splitter = re.split(r"[|\)\(>]+",e)
		print (splitter)
		prob = splitter[-1]
		if float(prob)<options.prune:
			continue
		anccod = splitter[0]
		ancaa = splitter[1]
		deraa = splitter[2]
		sub = {"anccod":anccod,"ancaa":ancaa,"deraa":deraa,"prob":prob,"text":e}
		if nname not in events:
			events[nname] = []
		events[nname].append(sub)
		if ancaa in counts:
			counts[ancaa] = counts[ancaa]+1
		else:
			counts[ancaa] = 1
		if deraa in counts:
			counts[deraa] = counts[deraa]+1
		else:
			counts[deraa] = 1
	
colors = {}
i = 0
for der, c in sorted(counts.items(),  key=lambda kv: (-kv[1], kv[0])):
		colors[der] = Paired_12.hex_colors[i]
		print(der)
		print(c)
		print(colors[der])
		i = i+1

#for der in set(derived.values()):
#		colors[der] = Paired_12.hex_colors[i]
#		print(der)
#		print(colors[der])
#		i = i+1


		
## todo : traverse and color grey
		
# for e in subtrees:		
	# if not options.large:
		# eventnodes = []
		# allnodes = []
		# for n in t.traverse():
		  # if not n.is_root():
				  # allnodes.append(n)  
		  # if n.name in e["events"]:
				  # eventnodes.append(n)
		# for n in allnodes:
		  # nname = n.name
		  # newchild = splitNode(n)
		  # if nname in e["subtree"] and not nname in e["events"]:
				  # e["subtree"].append(newchild.name) 




# Creates an independent node style for each node, which is
# initialized with a red foreground color.
t.ladderize()

#treecolor = "lightseagreen" #"mediumseagreen" # "lightseagreen" # "mediumaquamarine" #"seagreen" "lightseagreen"
treecolor = "lightgray" # "silver" #"DarkGray"
#eventscolor = "#006060" #"#707070" # "#006060"  #"teal" "Black"
for n in t.traverse():
	nstyle = NodeStyle()
	nstyle["hz_line_color"] = treecolor
	nstyle["vt_line_color"] = treecolor
	nstyle["size"] = 0
	nstyle["vt_line_width"] = 1
	nstyle["hz_line_width"] = 1
	n.set_style(nstyle)
	if n.name in events:
		for event in events[n.name]:
			anccolor = colors[event["ancaa"]]
			dercolor = colors[event["deraa"]]
			face = TextFace(event["text"],fgcolor = dercolor)
			face.background.color = anccolor
			n.add_face(face, column=1)

for aa in colors.keys():
	ts.legend.add_face(CircleFace(10, colors[aa]), column=0)
	ts.legend.add_face(TextFace(aa), column=1)
	
out=options.output
if options.prune:
	out=out + "_" + str(options.prune)
t.render(out + ".png", tree_style=ts, dpi=300, w=options.width, units="mm")
t.render(out + ".svg", tree_style=ts, dpi=300, w=options.width, units="mm")
