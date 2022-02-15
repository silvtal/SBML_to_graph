#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 12:12:15 2021

@author: silvtal
"""
from reframed import load_cbmodel
from reframed import Environment, CAFBA
from carveme.reconstruction.utils import load_media_db
import copy
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

def create_spare_medium(f, medium, media_db, init_env=None, outputdir=".", outputname="tempmedium.tsv", include_original=True, to_remove=None, objective="Growth"):
    """ Spare medium creation. This function takes a given SBML model and simulates its growth on a given medium. Then returns a list of metabolites that are the "spent media", that is, the original media plus the model's secretions and without the nutrients that the model has fully consumed.
    
    Parameters
    ----------
    f : str
        model file name
    medium : str
             medium name in the media_db
    media_db : str
               .tsv file with the medium composition
    init_env : Environment, optional
               ReFramed-generated Environment object (default: None)
    outputdir : str, optional
                output directory for the spent media file (default: current working directory)
    include_original : bool, optional
                       whether to include original medium compounds on the output file, except those specified in to_remove (default: None)
    to_remove : list, optional
                list of compounds to remove from the original medium (default: None)
    objective : str, optional
                name of objective function (default: "Growth")
    """
    # load and solve the spent_model
    s_model = load_cbmodel(f,flavor="fbc2")   
    if init_env==None:
        init_env = Environment.from_compounds(media_db[medium])
    
    init_env.apply(s_model)
    solution = CAFBA(s_model)  
    
    # select all positive flux reactions
    df=solution.to_dataframe()
    positive_fluxes= df[df["value"]>0]
    
    # pick all exchange reactions except sinks
    my_r= [r for r in s_model.get_exchange_reactions() if "sink" not in r]
    
    spare=[r for r in list(positive_fluxes.index) if r in my_r]
    
    # get metabolites from those reactions
    my_m=["".join(r.split("_")[2:-1]) for r in spare]

    # create file
    if include_original:
        if to_remove==None:
            my_m= media_db[medium]+my_m
        else:
            medcopy=copy.copy(media_db[medium])
            for c in to_remove:
                medcopy.remove(c)
                my_m= medcopy+my_m
    
    temp_media=pd.DataFrame({
        "medium" : ["temp"]*len(my_m),
        "description" : ["spent medium without original"]*len(my_m),
        "compound" : my_m
        })
    
    temp_media.to_csv( outputdir+"/"+outputname,sep="\t")



def metabolite_translator(names_table):
    """
    Metabolite name translator
    
    Parameters
    ----------
    names_table : str
                  text file with BiGG ID-metabolite
    
    Returns
    -------
    dict
         keys are bigg IDs with one metabolite name as unique value
    """

    names_table = pd.read_csv(names_table, sep="\t")
    code_to_name = dict()
    for index, row in names_table.iterrows():
        pretty_name = str(row['name']).replace(" ","_")
        compartment = row['bigg_id'].split("_")[-1]
        code_to_name[row['bigg_id']]=pretty_name+"_"+compartment
                          
    return(code_to_name)
    

    

def reaction_parser(all_reac={}, weights={},code_to_name=None,
                    exclude_reversible=False, 
                    exclude_irreversible=False, min_flux=0, 
                    given_color="black"):
    adj_list=[]
    for reac in range(len(all_reac)): 
        r=str(all_reac[reac]).split(":")[1]
        rname=str(all_reac[reac]).split(":")[0]
        flux=weights[list(weights.keys())[reac]]
        if abs(flux)<min_flux:
            continue # skip this iteration/reaction
        # parseo (diferente según reversibilidad)
        if "<->" in r and not exclude_reversible:
            s_and_p = r.split("<->")
            # decido cómo pintar su flecha:
            if flux>0: 
                ar="-|>"
            elif flux<0:
                ar="<|-"
            else:
                print("The following reaction's weight is zero...? Something must be wrong with your filtering")
                print(reac)
                continue
        elif "-->" in r and not exclude_irreversible:
            ar="-|>"
            s_and_p = r.split("-->")
        else:
            continue #skip excluded reactions
            
        if len(s_and_p)>1: # solo si es una reacción con sustrato y producto
            try:    
                sust=[s.replace(" ","").split("M_",maxsplit=1)[1] for s in s_and_p[0].split("+")]
                prod=[p.replace(" ","").split("M_",maxsplit=1)[1] for p in s_and_p[1].split("+")]
             #   print("sust",sust)
             #   print("prod",prod)
                print(rname+":"+str(sust)+"-->"+str(prod))
                for j in range(len(prod)):
                    for i in range(len(sust)):
                        new_edge=str(code_to_name[sust[i]]+" "+code_to_name[prod[j]]+ #weight: valor absoluto flujo
                       " {'weight': "+str(flux)+
                        
                       ", 'given_color': '"+given_color+"'"+
                       ", 'name':'"+ar+"'}")
                        adj_list.append(new_edge)
#                        print(code_to_name[sust[i]])
#                        print(ar)
#                        print(code_to_name[prod[j]])
#                        print("\n")

            except:
                print("Didn't know how to parse this: "+r)
                
              #  print(list(weights.keys())[reac])
              #  print(weights[list(weights.keys())[reac]])

    return(adj_list)



def my_draw(my_model, my_reactions, code_to_name, 
            exi, exr, minf, my_figsize,
            fntsize, textrot,
            adj_list=None, 
            given_color="black", 
            given_node_size=20,
            good_pos=None,
            label_options=None,
            highlight=None):    
    """
    Creates a graph. Calls reaction_parser.

    
    Parameters (REQUIRED)
    ----------
    my_model : CBModel 
               loaded with reframed.load_cbmodel
    my_reactions : dict
                   reactions and their weights. It comes directly from ReFramed's solution
                   example:
                      a = reframed_solution
                      my_reactions = {v:a.values[v] for v in a.values.keys()}
    code_to_name : dict 
                   created by metabolite_translator
    exi : bool
          whether to exclude or not irreversible functions from the graph
    exr : bool
          whether to exclude or not reversible functions from the graph
    minf : float
           minimum flux of reactions to be included in the graph
    my_figsize : tuple 
                 graph sizes (height and width)
    fntsize : int
              Font size for the labels
    textrot : float
              text rotation for the labels (degrees)
    
    Parameters (optional)
    ---------
    adj_list : list
               previously created graph. Useful for painting a graph on top of another
    given_color : str
                  color for edges
    given_node_size : int
                      size for nodes
    good_pos : dict
               node positions. Can be taken from another graph with 
               networkx.drawing.nx_agraphgraphviz_layout(existing_graph, prog='neato')
    label_options : dict
                    example: "{"ec": "black", "fc": "white", "alpha": 0.7}", where "ec"
                    means box border color and "fc" means box fill color
    highlight :  dict
                 the dictionary keys are a color and the corresponding value is a list
                 of metabolites to highlight
               
    """
    
    # Create a graph with reaction_parser and code_to_name
    ###############
    all_reac=[]
    for r in my_reactions.keys():
        all_reac.append(my_model.reactions[r])
    if adj_list == None:
        adj_list = reaction_parser(all_reac, weights=my_reactions, code_to_name=code_to_name,
                                 exclude_irreversible=exi, exclude_reversible=exr, min_flux=minf, 
                                 given_color=given_color)
    else:
        adj_list += reaction_parser(all_reac, weights=my_reactions, code_to_name=code_to_name,
                                 exclude_irreversible=exi, exclude_reversible=exr, min_flux=minf, 
                                 given_color=given_color)
    G = nx.parse_edgelist(adj_list, nodetype = str,data=True, create_using=nx.DiGraph)


    # layouts
    ##########
    G_4_layout = nx.parse_edgelist(adj_list, nodetype = str, data = True) 

    df = pd.DataFrame(index=G_4_layout.nodes(), columns=G_4_layout.nodes())
    for row, data in nx.shortest_path_length(G_4_layout):
        for col, dist in data.items():
            df.loc[row,col] = dist

    df = df.fillna(df.max().max())
     
    if good_pos == None:
        pos=nx.kamada_kawai_layout(G, dist=df.to_dict())
        label_pos=nx.kamada_kawai_layout(G, dist=df.to_dict())
    else:
        pos = good_pos
        label_pos = good_pos

    # draw
    # ====
    fig=plt.figure(figsize=my_figsize, dpi= 200, facecolor='w', edgecolor='k')

    nx.draw_networkx_nodes(G, pos,
                           node_color="darkgrey",
                           node_size=given_node_size)    

    ###label_options = {"ec": "k", "fc": "white", "alpha": 0.7}
    text=nx.draw_networkx_labels(G, label_pos, font_size=fntsize, clip_on=False, bbox=label_options)
    ### optional highlight:
    if (highlight != None):
        if (label_options == None):
            label_options = {}
        for k in highlight.keys():
            selected_nodes =  [n for n,v in G.nodes(data="True") if n in highlight[k]]
            label_options["ec"] = k
            label_options["fc"] = k
            label_options["alpha"] = 0.7
            text=nx.draw_networkx_labels(G.subgraph(selected_nodes), label_pos, font_size=fntsize, clip_on=False, bbox=label_options)

    for _,t in text.items():
        t.set_rotation(textrot)

    for edge in G.edges(data=True):
        ar = edge[2]['name']
        nx.draw_networkx_edges(G, pos,
                               node_size=given_node_size,
                               edgelist=[(edge[0],edge[1])],
                               arrowstyle=str(ar),
                               arrowsize=14,
                               edge_color=edge[2]["given_color"],
                               width=1)

    ax = plt.gca()
    ax.set_axis_off()
    return(ax, adj_list)
