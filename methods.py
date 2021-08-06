import semopy
from semopy import Model
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
from matplotlib import pyplot as plt
#read and process file
pathways=pd.read_csv("data/pathways.tsv", sep='\t')
esets=pd.read_csv("data/controls_counts_norm.csv")
esets.index=esets["Unnamed: 0"].values
esets=esets.drop(columns="Unnamed: 0")
gene_edges=pd.read_csv("data/gene_edges.tsv", sep='\t')



#given a selected pathway, return the processed edges of the pathway
def read_pathway(pathway_name):
    pathway_id=np.where(pathways["pathway_name"]==pathway_name)[0]
    apathway=pathways["nodes"].loc[int(pathway_id)]
    apathway=apathway.split(";")
    apathway=[int(x) for x in apathway]
    pathway_edges=pd.DataFrame([list(x) for x in gene_edges.values if x[0] in apathway])
    pathway_edges=pd.DataFrame([x for x in pathway_edges.values if str(x[0]) in esets.index])
    pathway_edges=pd.DataFrame([x for x in pathway_edges.values if str(x[1]) in esets.index])
    return pathway_edges

#calculates the adj_matrix and the dictionary to find the nodes, with new names, in the matrix
def build_adj(pathway_edges):
    pathway_edges_0=pathway_edges[0].unique()
    pathway_edges_1=pathway_edges[1].unique()
    nodes=list(np.hstack((pathway_edges_0,pathway_edges_1)))
    nodes_renamed={}
    inv_nodes_renamed={}
    for e,x in enumerate(nodes):
        nodes_renamed[x]=e
        inv_nodes_renamed[e]=x
    nodes=len(nodes)
    adj_matrix=np.zeros((nodes,nodes))
    for x in pathway_edges.values:
        adj_matrix[nodes_renamed[x[0]]][nodes_renamed[x[1]]]=x[2]
    return (adj_matrix,nodes_renamed,inv_nodes_renamed)

#gets the triads
def get_triad(G):
    all_cliques= nx.enumerate_all_cliques(G)
    triad_cliques=[x for x in all_cliques if len(x)==3 ]
    return triad_cliques


#does the SEM analysis, weight every edge and remove any invalid triad in the list of cliques
def calculate_weighted_edges(triad_cliques, adj_matrix,inv_nodes_renamed):
    weighted_edges={}
    new_triad_cliques=[]
    first_label=""
    second_label=""
    third_label=""
    #SEM model
    mod = """ 
              y ~ x1 + x2
              """
    #for each triad
    for triad in range(len(triad_cliques)):
        triad_matrix=np.zeros((3,3))
        for i,x in enumerate(triad_cliques[triad]):
            for j,y in enumerate(triad_cliques[triad]):
                triad_matrix[i][j]=adj_matrix[x][y]
        zeros_count=np.array([len(np.where(x==0)[0]) for i,x in enumerate(triad_matrix) ])
        #filtering any invalid clique
        if (sum(zeros_count)==6):
            new_triad_cliques.append(triad_cliques[triad])
        else:
            continue
        #the first is always the rows of the matrix with two edges (one 0 in the matrix) to the other genes
        #the same goes for the other two genes (one edge for the second gene and zero for the third)
        first_index=int(np.where(zeros_count==1)[0])
        second_index=int(np.where(zeros_count==2)[0])
        third_index=int(np.where(zeros_count==3)[0])
        #gets the gene expression of each gene
        first_label=str(inv_nodes_renamed[triad_cliques[triad][first_index]])
        second_label=str(inv_nodes_renamed[triad_cliques[triad][second_index]])
        third_label=str(inv_nodes_renamed[triad_cliques[triad][third_index]])
        first_gene=(list(esets.loc[first_label,:].values),0)
        second_gene=(list(esets.loc[second_label,:].values),1)
        third_gene=(list(esets.loc[third_label,:].values),2)

        y=third_gene
        x1=first_gene
        x2=second_gene
        to_df={"y":y[0],"x1":x1[0],"x2":x2[0]}
        data=pd.DataFrame(to_df).replace(np.inf, np.nan).replace(-np.inf, np.nan).dropna()
        m = Model(mod)
        #SEM analysis
        r = m.fit(data)
        #sum of latent variables found by SEM
        fac_sum=np.abs(r.x[0]+r.x[1])
        #check if the edge between the first and the third node is less significant than the others two, this edge will be weighted as 1,
        #the other two as 0
        if (np.abs(r.x[0])<fac_sum*0.1):
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],0))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],0)]
            if (second_label+","+third_label in weighted_edges) :
                weighted_edges[second_label+","+third_label].append((r.x[1],1))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],1)]

            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],1))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],1)]
        #the other two are less significant
        elif(np.abs(r.x[1])<fac_sum*0.1):
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],1))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],1)]
            if (second_label+","+third_label in weighted_edges):  
                weighted_edges[second_label+","+third_label].append((r.x[1],0))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],0)]
            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],0))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],0)]
        #there are no most significant edge in the clique, all edges will be weighted as -1
        else:
            if (first_label+","+third_label in weighted_edges):
                weighted_edges[first_label+","+third_label].append((r.x[0],-1))

            else:
                weighted_edges[first_label+","+third_label]=[(r.x[0],-1)]
            if (second_label+","+third_label in weighted_edges):  
                weighted_edges[second_label+","+third_label].append((r.x[1],-1))
            else:
                weighted_edges[second_label+","+third_label]=[(r.x[1],-1)]
            if(first_label+","+second_label in weighted_edges) :
                weighted_edges[first_label+","+second_label].append((r.x[1],-1))
            else:
                weighted_edges[first_label+","+second_label]=[(r.x[1],-1)]
    return weighted_edges, new_triad_cliques

def evaluate_edges(weighted_edges):
    to_remove=[]
    equi_values={}
    essential_edges=[]
    #for each weighted edge
    for x in weighted_edges.items():
        zeros=0
        ones=0
        minus=0
        for z in x[1]:
            if (z[1]==0):
                zeros+=1
            elif (z[1]==1):
                ones+=1
            else:
                minus+=1
        if (ones==0):
            if (minus==0):
                to_remove.append(x[0])
            else:
                m=(zeros+minus)/2
                #remove crition
                if ((minus+zeros)/(zeros*minus+1)*zeros/(minus+1)>((m+m)/(m*m+1))*m/(m+1)):
                    to_remove.append(x[0])
        else:
            essential_edges.append(x[0])
        if (ones==0):
            equi_values[x[0]]=round((minus+zeros)/(zeros*minus+1)*(zeros)/(minus+1),3)
        else:
            equi_values[x[0]]=0
    return to_remove, equi_values, essential_edges

