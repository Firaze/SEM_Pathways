import streamlit as st
import streamlit.components.v1 as components
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network
from functions import *
import PIL
#Network(notebook=True)

# makes Network show itself with repr_html
pathways_name=pd.read_csv("data/pathways.tsv", sep='\t')["pathway_name"]

#excange of pathway to make easier the findings of the (interesting) pathway on the 25th position 
tmp=pathways_name.iloc[0]
pathways_name.iloc[0]=pathways_name.iloc[25]
pathways_name.iloc[25]=tmp

st.set_page_config(layout="wide")

#populating sidebar
st.sidebar.title('Choose a pathway')
option=st.sidebar.selectbox('',pathways_name)
normal_edges=st.sidebar.checkbox('Show expression/suppression edges')
removed_edges=st.sidebar.checkbox('Show removed edges',value=True)
st.sidebar.text("Edge legend:")
w=25
h=5
red = np.zeros((h, w, 3), dtype=np.uint8)
red[0:]=[255, 0, 0]
green = np.zeros((h, w, 3), dtype=np.uint8)
green[0:]=[0, 128, 0]
blue = np.zeros((h, w, 3), dtype=np.uint8)
blue[0:]=[30,144,255]
yellow = np.zeros((h, w, 3), dtype=np.uint8)
yellow[0:]=[255, 255,0]
st.sidebar.image(blue, caption='Expression edges')
st.sidebar.image(yellow, caption='Suppression edges')
st.sidebar.image(green, caption='Part of triad edges')
st.sidebar.image(red, caption='Removed edges')


#skip_calcs if an invalid pathways was selected
skip_calcs=False
#selects pathway
pathway_edges=read_pathway(option)
#checks if the pathway is invalid
if (len(pathway_edges)==0):
    #option=last_selection
    skip_calcs=True
    st.error("Edges not found, try another pathway!")
else:
    skip_calcs=False
    #last_selection=option
    pathway_edges=read_pathway(option)
    #calculates the adj_matrix and the dictionary to find the nodes, with new names, in the matrix
    adj_matrix,nodes_renamed,inv_nodes_renamed=build_adj(pathway_edges)
    G = nx.from_numpy_matrix(adj_matrix)
    #gets the triads
    triad_cliques=get_triad(G)
    #does the SEM analysis, weight every edge and remove any invalid triad in the list of cliques
    weighted_edges,triad_cliques=calculate_weighted_edges(triad_cliques, adj_matrix,inv_nodes_renamed)
    #calculate which edges must be removed, the edges that must be kept and all their significative values 
    to_remove, equi_values, essential_edges=evaluate_edges(weighted_edges)
    relabel={}
    #relabel the nodes of the graph to obtain the real names of the genes
    for e,node in enumerate( G.nodes()):
        relabel[e]=str(inv_nodes_renamed[node])
    net=Network(height="825px",notebook=True,directed=True,width="1800px", bgcolor='#222222', font_color='white')
    
    #filters the nodes if the user wants to show only triad results
    triad_nodes=set()
    _=[triad_nodes.add(str(inv_nodes_renamed[y])) for x in triad_cliques for y in x]
    triad_nodes=list(triad_nodes)
    
    #adds nodes in the grapg
    for i,node in relabel.items():
        if normal_edges:
            net.add_node(str(node))
        #add only triad nodes if the user wants triad related results
        elif node in triad_nodes:
            net.add_node(str(node))
    #same for the edges
    if (normal_edges):
        for edge in pathway_edges.values:
                if(edge[2]==-1):
                    net.add_edge(str(edge[0]), str(edge[1]), color="yellow")
                else:
                    net.add_edge(str(edge[0]), str(edge[1]))
    #colors the edge of green if it is not removed, red otherwise
    #also assigns the values, calculated some lines above, to each edge
    for triad in triad_cliques:
        for i,x in enumerate(triad):
            for j,y in enumerate(triad):
                value=""
                isessential=""
                tmp=pathway_edges[(pathway_edges[0]==inv_nodes_renamed[triad[i]]) & (pathway_edges[1]==inv_nodes_renamed[triad[j]])].values
                if (len(tmp)>0):
                    start_node,to_node,weight=tmp[0]
                else:
                    continue
                if ((str(start_node)+","+str(to_node)) in to_remove): 
                    if (removed_edges==False):
                        continue
                    color="red"
                    size=3
                    value+=", equilibrium factor:  "+str(equi_values[str(start_node)+","+str(to_node)])
                else:
                    color="green"
                    size=3
                    value+=", equilibrium factor:  "+str(equi_values[str(start_node)+","+str(to_node)])
                if ((str(start_node)+","+str(to_node)) in essential_edges):   
                    isessential="Essential "
                if (weight==1):
                    net.add_edge(str(start_node), str(to_node), color=color, width=size,title=isessential+"Expression edge"+value)
                else:
                    net.add_edge(str(start_node), str(to_node), color=color, width=size,title=isessential+"Suppression edge"+value)
    #sets the garph
    net.hrepulsion(node_distance=120, central_gravity=0.0, spring_length=100, spring_strength=0, damping=0.09)
    net.show("data/graph.html")
    #reads and shows the graph
    HtmlFile = open("data/graph.html", 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code, height = 850,width=1850)

