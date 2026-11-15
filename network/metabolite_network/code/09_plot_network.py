from pathlib import Path
from plot_PPI import plot_ppi
from plot_KYNU import plot_kynu
from figure_style import COLORS, save
import textwrap
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.text import Text

def curve(ax,x0,y0,x1,y1,color,width=1,alpha=.5):
    mid=(x0+x1)/2
    path=MplPath([(x0,y0),(mid,y0),(mid,y1),(x1,y1)],[MplPath.MOVETO,MplPath.CURVE4,MplPath.CURVE4,MplPath.CURVE4])
    ax.add_patch(PathPatch(path,fill=False,color=color,lw=width,alpha=alpha))

def main(root):
    result=root/'results'
    sources=root/'data/figure_source'
    sources.mkdir(exist_ok=True)
    nodes=pd.read_csv(result/'network_node_statistics.csv').set_index('gene_symbol')
    edges=pd.read_csv(result/'STRING_PPI_edges.csv')
    graph=nx.from_pandas_edgelist(edges,'gene_A','gene_B',edge_attr='score')
    plot_ppi(root,nodes,graph,COLORS,save)
    plot_kynu(root)
    modules=pd.read_csv(result/'MCODE_modules.csv')
    fig=plt.figure(figsize=(7,7.2))
    grids=fig.add_gridspec(2,2,left=.06,right=.97,top=.94,bottom=.09,hspace=.27,wspace=.25)
    for i in range(3):
        ax=fig.add_subplot(grids[i//2,i%2]); row=modules.iloc[i]
        genes=row.genes.split(';'); sub=graph.subgraph(genes)
        pos=nx.circular_layout(sorted(sub))
        nx.draw_networkx_edges(sub,pos,ax=ax,width=.5,edge_color='#ABB3BA',alpha=.40)
        nx.draw_networkx_nodes(sub,pos,ax=ax,node_size=65,node_color=[COLORS['pink'] if nodes.loc[g,'consensus_hub'] else COLORS['purple'] for g in sub],edgecolors='none')
        for gene,(x,y) in pos.items():
            angle=np.degrees(np.arctan2(y,x))
            rotation=angle if x>=0 else angle+180
            ax.text(x*1.12,y*1.12,gene,fontsize=7,ha='left' if x>=0 else 'right',va='center',rotation=rotation,rotation_mode='anchor')
        ax.set_xlim(-1.75,1.75);ax.set_ylim(-1.9,1.95);ax.set_aspect('equal');ax.axis('off')
        ax.set_title(f'MCODE module {row.module}\n{row["size"]} genes; score = {row.score:.2f}',pad=4)
        ax.text(-.05,1.09,chr(97+i),transform=ax.transAxes,weight='bold',fontsize=11)
    ax=fig.add_subplot(grids[1,1]); top=nodes.sort_values(['MCC'],ascending=False).head(10).iloc[::-1]
    ax.barh(top.index,np.log10(top.MCC+1),color=[COLORS['pink'] if x else COLORS['purple'] for x in top.consensus_hub],height=.65)
    ax.set_xlabel('log$_{10}$(MCC + 1)');ax.set_title('cytoHubba top 10');ax.text(-.18,1.09,'d',transform=ax.transAxes,weight='bold',fontsize=11)
    top.to_csv(sources/'MCC_display.csv')
    save(fig,root,'Fig6_2_modules_and_hubs')
    enrich=pd.read_csv(result/'Metascape_enrichment_FDR05.csv')
    fig,axes=plt.subplots(3,1,figsize=(7,8.6),gridspec_kw={'height_ratios':[1,1,1.35]})
    fig.subplots_adjust(left=.60,right=.96,top=.94,bottom=.09,hspace=.60)
    categories=[('GO Biological Processes','Biological process',COLORS['teal']),('GO Cellular Components','Cellular component',COLORS['orange']),('GO Molecular Functions','Molecular function',COLORS['purple'])]
    displayed=[]
    for i,(category,title,color) in enumerate(categories):
        ax=axes[i]; subset=enrich[enrich.Category.eq(category)].sort_values(['q_value','GO']).head(5).iloc[::-1]
        displayed.append(subset)
        labels=['\n'.join(textwrap.wrap(x,43)) for x in subset.Description]
        ax.barh(np.arange(len(subset)),-subset.log10_q,color=color,height=.65)
        ax.set_yticks(np.arange(len(subset)),labels)
        ax.set_title(title,loc='left',pad=8)
        ax.set_xlabel('-log$_{10}$(BH-adjusted P)')
        ax.text(-1.55,1.13,chr(97+i),transform=ax.transAxes,weight='bold',fontsize=11)
    pd.concat(displayed).to_csv(sources/'GO_display_terms.csv',index=False)
    save(fig,root,'Fig6_3_GO_enrichment')
    links=pd.read_csv(result/'KEGG_display_gene_pathway_edges.csv')
    terms=links[['pathway_id','pathway_name','q_value']].drop_duplicates().reset_index(drop=True)
    genes=sorted(links.gene_symbol.unique(),key=lambda g:(links.loc[links.gene_symbol.eq(g),'pathway_id'].iloc[0],g))
    gy=dict(zip(genes,np.linspace(.96,.04,len(genes))))
    fig,ax=plt.subplots(figsize=(7,9));fig.subplots_adjust(left=.03,right=.97,top=.96,bottom=.03)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    palette=[COLORS['teal'],COLORS['pink'],COLORS['purple']]
    for i,row in terms.iterrows():
        assigned=links.loc[links.pathway_id.eq(row.pathway_id),'gene_symbol'].tolist()
        center=.82-i*.32
        endpoints=np.linspace(center+.08,center-.08,len(assigned))
        ax.add_patch(Rectangle((.69,center-.085),.018,.17,color=palette[i],ec='none'))
        label='\n'.join(textwrap.wrap(row.pathway_name,20))
        ax.text(.73,center,label+f'\nq = {row.q_value:.2g}',va='center',fontsize=9,linespacing=1.6)
        for gene,y1 in zip(assigned,endpoints):
            curve(ax,.32,gy[gene],.69,y1,palette[i],2.3,.42)
    for gene,y in gy.items():
        ax.text(.285,y,gene,ha='right',va='center',fontsize=8)
        ax.add_patch(Rectangle((.305,y-.004),.015,.008,facecolor=COLORS['ink'],edgecolor='none'))
    ax.text(.31,.995,'Mapped genes',ha='right',weight='bold',fontsize=10)
    ax.text(.69,.995,'Top KEGG pathways',weight='bold',fontsize=10)
    links.to_csv(sources/'KEGG_Sankey_edges.csv',index=False)
    save(fig,root,'Fig6_4_KEGG_gene_pathway_Sankey')
    bridge=pd.read_csv(result/'KP_expression_bridge_mapping.csv')
    evidence=pd.read_csv(result/'KP_expression_evidence_matrix.csv')
    gene_order=['KYNU','SLC7A5','SLC7A8','ALAS1','MBOAT1','GPAT3','ABCG2']
    compound_order=['Alanine','5-Aminolevulinic acid','1-Octadecanoyl-sn-glycero-3-phosphoethanolamine','1-Palmitoyllysophosphatidate','Dehydroepiandrosterone sulfate','Indoxyl sulfate','testosterone sulfate']
    bridge['metabolite_name']=bridge.metabolite_name.str.strip()
    genes_y=dict(zip(gene_order,np.linspace(.85,.15,7)))
    compounds_y=dict(zip(compound_order,[.85,.63,.49,.36,.23,.12,.03]))
    datasets=['GSE272548_CAS_90-15-3','GSE272548_CAS_135-19-3','GSE274318_CAS_90-15-3','GSE274318_CAS_135-19-3']
    fig,ax=plt.subplots(figsize=(7,5.3));fig.subplots_adjust(left=.02,right=.98,top=.94,bottom=.15)
    ax.set_xlim(0,1.02);ax.set_ylim(-.05,1.1);ax.axis('off')
    for row in bridge[['metabolite_name','gene_symbol']].drop_duplicates().itertuples():
        curve(ax,.34,compounds_y[row.metabolite_name],.53,genes_y[row.gene_symbol],COLORS['orange'] if row.gene_symbol=='KYNU' else COLORS['grey'],1.7,.9)
    for name,y in compounds_y.items():
        ax.text(.32,y,'\n'.join(textwrap.wrap(name,30)),ha='right',va='center',fontsize=8,linespacing=1.1)
        ax.scatter([.34],[y],s=42,color=COLORS['purple'],edgecolors='none',zorder=3)
    for gene,y in genes_y.items():
        ax.scatter([.53],[y],s=120 if gene=='KYNU' else 65,color=COLORS['orange'] if gene=='KYNU' else COLORS['teal'],edgecolors='none',zorder=3)
        ax.text(.55,y,gene,ha='left',va='center',fontsize=8,weight='bold' if gene=='KYNU' else 'normal')
        for dataset,x in zip(datasets,[.73,.81,.90,.98]):
            record=evidence[(evidence.gene_symbol==gene)&(evidence.dataset_id==dataset)].iloc[0]
            if record.status=='responsive':
                ax.scatter([x],[y],s=70,color=COLORS['pink'],edgecolors='none')
            elif record.status=='not_responsive':
                ax.scatter([x],[y],s=60,facecolor='white',edgecolor='#8E9CA5',linewidth=.8)
            else:
                ax.scatter([x],[y],s=30,color='#8E9CA5',marker='x',linewidth=.8)
    ax.text(.17,1.08,'Differential metabolite annotations',ha='center',weight='bold')
    ax.text(.55,1.08,'Mapped genes',ha='center',weight='bold')
    ax.text(.77,1.08,'MCF7',ha='center',weight='bold');ax.text(.94,1.08,'U2OS',ha='center',weight='bold')
    for x,label in zip([.73,.81,.90,.98],['1-Nap','2-Nap','1-Nap','2-Nap']):
        ax.text(x,1,label,ha='center',fontsize=8,rotation=45)
    handles=[Line2D([],[],marker='o',ls='',color=COLORS['pink'],label='Exposure-responsive',markeredgewidth=0),Line2D([],[],marker='o',ls='',markerfacecolor='white',markeredgecolor='#8E9CA5',label='Not responsive'),Line2D([],[],marker='x',ls='',color='#8E9CA5',label='Not testable')]
    fig.legend(handles=handles,loc='lower center',ncol=3,bbox_to_anchor=(.5,.02),fontsize=8)
    bridge.to_csv(sources/'KP_expression_bridge_edges.csv',index=False)
    evidence.to_csv(sources/'KP_expression_matrix.csv',index=False)
    save(fig,root,'Fig6_5_KP_expression_bridge')
    print('Exported six figures in PDF, SVG, TIFF and PNG')

if __name__=='__main__':
    main(Path(__file__).resolve().parents[1])
