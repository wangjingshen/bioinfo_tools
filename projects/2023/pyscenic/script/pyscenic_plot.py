import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import umap
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
import natsort
from collections import OrderedDict
from math import ceil
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from adjustText import adjust_text

from pyscenic.plotting import plot_binarization
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss


class Pyscenic_plot:
    def __init__(self, args):
        self.pyscenic_loom = args.pyscenic_loom
        self.h5ad = args.h5ad
        self.threads = int(args.threads)
        self.ncol = int(args.ncol)
        self.outdir = args.outdir
        if not os.path.exists(self.outdir):
            os.system(f'mkdir -p {self.outdir}')
        #sc.settings.figdir(self.outdir)
    
    def add_metadata(self):
        # auc reduction --
        lf = lp.connect(self.pyscenic_loom, mode='r+', validate=False )
        auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
        # UMAP
        runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
        dr_umap = runUmap( auc_mtx )
        pd.DataFrame(dr_umap, columns=['X', 'Y'], index = auc_mtx.index).to_csv( f"{self.outdir}/pyscenic_auc_umap.txt", sep='\t')
        # tSNE
        tsne = TSNE( n_jobs = self.threads )
        dr_tsne = tsne.fit_transform( auc_mtx )
        pd.DataFrame(dr_tsne, columns=['X', 'Y'], index = auc_mtx.index).to_csv( f"{self.outdir}/pyscenic_auc_tsne.txt", sep='\t')


        # scenic output
        adata = sc.read_h5ad(self.h5ad)
        meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
        #exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
        #auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
        dr_umap = pd.read_csv( f'{self.outdir}/pyscenic_auc_umap.txt', sep='\t', header=0, index_col=0 )
        dr_tsne = pd.read_csv( f'{self.outdir}/pyscenic_auc_tsne.txt', sep='\t', header=0, index_col=0 )
        auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
        regulons = lf.ra.Regulons
        regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
        # regulon thresholds
        rt = meta['regulonThresholds']
        for i,x in enumerate(rt):
            tmp = x.get('regulon').replace("(","_(")
            x.update( {'regulon': tmp} )
    

        tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])
        Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
        Embeddings_X = pd.concat( [
            pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
            pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
            dr_tsne['X'] ,
            dr_umap['X']
        ], sort=False, axis=1, join='outer' )
        Embeddings_X.columns = ['1','2','3','4']

        Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
        Embeddings_Y = pd.concat( [
            pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
            pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
            dr_tsne['Y'] ,
            dr_umap['Y']
        ], sort=False, axis=1, join='outer' )
        Embeddings_Y.columns = ['1','2','3','4']

        ### metadata
        metaJson = {}

        metaJson['embeddings'] = [
        {"id": -1, "name": f"Scanpy t-SNE (highly variable genes)"},
        {"id": 1, "name": f"Scanpy UMAP  (highly variable genes)"},
        {"id": 2, "name": "Scanpy PC1/PC2"},
        {"id": 3, "name": "SCENIC AUC t-SNE"},
        {"id": 4, "name": "SCENIC AUC UMAP"}]

        metaJson["clusterings"] = [
            {"id": 0,
             "group": "Scanpy",
             "name": "Scanpy default resolution",   # Scanpy louvain default resolution
             "clusters": []}
        ]

        metaJson["metrics"] = [
            {"name": "nUMI"}, 
            {"name": "nGene"}, 
            {"name": "Percent_mito"}
        ]

        metaJson["annotations"] = [
            { "name": "clusters_Scanpy",
              "values": list(set( adata.obs['cluster'].astype(np.str) ))},
        ]

        # SCENIC regulon thresholds:
        metaJson["regulonThresholds"] = rt

        for i in range(max(set([int(x) for x in adata.obs['cluster']])) + 1):   # for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
            clustDict = {}
            clustDict['id'] = i
            clustDict['description'] = f'Unannotated Cluster {i + 1}'
            metaJson['clusterings'][0]['clusters'].append(clustDict)
    
        clusterings = pd.DataFrame()
        clusterings["0"] = adata.obs['cluster'].values.astype(np.int64)    # clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)

        def dfToNamedMatrix(df):
            arr_ip = [tuple(i) for i in df.values]
            dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
            arr = np.array(arr_ip, dtype=dtyp)
            return arr

        col_attrs = {
            "CellID": np.array(adata.obs.index),
            "cell_type": np.array(adata.obs['cell_type']),
            "nUMI": np.array(adata.obs['n_counts'].values),
            "nGene": np.array(adata.obs['n_genes'].values),
            #"Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
            #"Genotype": np.array(adata.obs['Genotype'].values),
            #"Timepoint": np.array(adata.obs['Timepoint'].values),
            #"Sample": np.array(adata.obs['Sample'].values),
            "Percent_mito": np.array(adata.obs['percent_mito'].values),
            "Embedding": dfToNamedMatrix(tsneDF),
            "Embeddings_X": dfToNamedMatrix(Embeddings_X),
            "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
            "RegulonsAUC": dfToNamedMatrix(auc_mtx),
            "Clusterings": dfToNamedMatrix(clusterings),
            "ClusterID": np.array(adata.obs['cluster'].values)
        }

        row_attrs = {
            "Gene": lf.ra.Gene,
            "Regulons": regulons,
        }

        attrs = {
            "title": "sampleTitle",
            "MetaData": json.dumps(metaJson),
            "Genome": 'hg38',
            "SCopeTreeL1": "",
            "SCopeTreeL2": "",
            "SCopeTreeL3": ""
        }

        # compress the metadata field:
        attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

        lp.create(filename = f'{self.outdir}/pyscenic_output_add_metadata.loom', layers = lf[:,:], row_attrs = row_attrs, col_attrs = col_attrs, file_attrs = attrs)
        lf.close() # close original pyscenic loom file


    def plot(self):
        # scenic output
        lf = lp.connect( f'{self.outdir}/pyscenic_output_add_metadata.loom', mode='r', validate=False )
        meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
        exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
        auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
        # create a dictionary of regulons:
        regulons = {}
        for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
            regulons[i] =  list(r[r==1].index.values)
    
        # cell annotations from the loom column attributes:
        cellAnnot = pd.concat(
            [
                pd.DataFrame( lf.ca.cell_type, index=lf.ca.CellID ),
                pd.DataFrame( lf.ca.ClusterID, index=lf.ca.CellID ),
                #pd.DataFrame( lf.ca.Louvain_clusters_Scanpy, index=lf.ca.CellID ),
                pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
                pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
                pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
            ],
            axis=1
        )
        cellAnnot.columns = [
            'cell_type',
            'ClusterID',
            #'Louvain_clusters_Scanpy',
            'Percent_mito',
            'nGene',
            'nUMI']

        # capture embeddings:
        dr = [
            pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
        ]
        dr_names = [
            meta['embeddings'][0]['name'].replace(" ","_")
        ]

        # add other embeddings
        drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
        dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )

        for i in range( len(drx.columns) ):
            dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
            dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )

        # rename columns:
        for i,x in enumerate( dr ):
            x.columns = ['X','Y']
    
        lf.close()

        # set color dict
        color_dict = dict(zip( sorted(list(set(cellAnnot['cell_type']))), sns.color_palette(palette='bright', n_colors=len(set(cellAnnot['cell_type'])) )))

        # 
        def colorMap( x, palette='bright' ):
            #n=len(set(x))
            #cpalette = sns.color_palette(palette,n_colors=n )
            #cdict = dict( zip( list(set(x)), cpalette ))
            cdict = color_dict
            cmap = [ cdict[i] for i in x ]
            cdict = OrderedDict( natsort.natsorted(cdict.items()) )
            return cmap,cdict

        def drplot( dr, colorlab, ax, palette='bright', title=None, **kwargs ):
            cmap,cdict = colorMap( colorlab, palette )
            for lab,col in cdict.items():  
                ix = colorlab.loc[colorlab==lab].index
                ax.scatter( dr['X'][ix], dr['Y'][ix], s=10, c=[col]*len(ix), alpha=0.7, label=lab, edgecolors='none')
            if( title is not None ):
                ax.set_title(title, fontsize='x-large');
            #
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

        plt.rcParams.update({'font.size':12})

        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(10,5), dpi=50 )

        drplot( dr[1], colorlab=cellAnnot['cell_type'], ax=ax1, palette='bright', s=2, title='Highly variable genes - UMAP' )
        drplot( dr[4], colorlab=cellAnnot['cell_type'], ax=ax2, palette='bright', s=2, title='pyscenic AUC - UMAP' )
        ax2.legend(loc='right', bbox_to_anchor=(1.5, 0.5), ncol=1, markerscale=2, fontsize='medium', frameon=False, title="cell_type")

        plt.tight_layout()
        plt.savefig(f"{self.outdir}/pyscenic_cell_type_auc.pdf", dpi=600, bbox_inches = "tight")


        rss_cellType = regulon_specificity_scores( auc_mtx, cellAnnot['cell_type'] )
        cats = sorted(list(set(cellAnnot['cell_type'])))
        nrow = ceil(len(cats)/self.ncol)
        fig = plt.figure(figsize=(4*self.ncol, 4*nrow))
        for c,num in zip(cats, range(1,len(cats)+1)):
            if(num > len(cats)):
                break
            x=rss_cellType.T[c]
            ax = fig.add_subplot(nrow, self.ncol, num)
            plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
            ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
            for t in ax.texts:
                t.set_fontsize(8)
            ax.set_ylabel('')
            ax.set_xlabel('')
            adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )
 
        fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='medium')
        fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='medium')
        plt.tight_layout()
        plt.rcParams.update({
            'figure.autolayout': True,
            'figure.titlesize': 'medium' ,
            'axes.labelsize': 'medium',
            'axes.titlesize':'medium',
            'xtick.labelsize':'medium',
            'ytick.labelsize':'medium'
        })
        plt.savefig(f"{self.outdir}/pyscenic_plot_cellType_RSS_top5.pdf", dpi=600, bbox_inches = "tight")


        # heatmap --
        topreg = []
        for i,c in enumerate(cats):
            print(i,c,list(rss_cellType.T[c].sort_values(ascending=False)[:5].index))  # for check
            topreg.extend(
                list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
            )
        topreg = list(set(topreg))
        # Generate a Z-score for each regulon to enable comparison between regulons
        auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
        for col in list(auc_mtx.columns):
            auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
        #auc_mtx_Z.sort_index(inplace=True)
            
        def palplot(pal, names, colors=None, size=1):
            n = len(pal)
            f, ax = plt.subplots(1, 1, figsize=(n * size, size))
            ax.imshow(np.arange(n).reshape(1, n),
                cmap=mpl.colors.ListedColormap(list(pal)),
                interpolation="nearest", aspect="auto")
            ax.set_xticks(np.arange(n) - .5)
            ax.set_yticks([-.5, .5])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            colors = n * ['k'] if colors is None else colors
            for idx, (name, color) in enumerate(zip(names, colors)):
                ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
            return f

        colors = sns.color_palette('bright',n_colors=len(cats) )
        #colorsd = dict( zip( cats, colors ))
        colorsd = color_dict
        colormap = [ colorsd[x] for x in cellAnnot['cell_type'] ]

        sns.set()
        sns.set(font_scale=0.8)
        fig = palplot( colors, cats, size=1.0)
        plt.savefig(f"{self.outdir}/pyscenic_plot_heatmap_legends.pdf", bbox_inches = "tight")   # dpi=600

        sns.set(font_scale=2)
        g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
            yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
            cmap="YlGnBu", figsize=(18,16) )
        g.cax.set_visible(True)
        g.ax_heatmap.set_ylabel('')
        g.ax_heatmap.set_xlabel('')
        plt.savefig(f"{self.outdir}/pyscenic_plot_heatmap_top5.pdf", dpi=600, bbox_inches = "tight")    

    def run(self):
        self.add_metadata()
        self.plot()


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--pyscenic_loom', help='pyscenic loom', required=True)
    parsers.add_argument('--h5ad', help='h5ad', required=True)
    parsers.add_argument('--threads', help='threads', required=True)
    parsers.add_argument('--ncol', help='ncol', required=True)
    parsers.add_argument('--outdir', help='outdir', required=True)
    args = parsers.parse_args()
    if(args.ncol == None):
        args.ncol = 3

    runner = Pyscenic_plot(args) 
    runner.run()

if __name__ == '__main__':
    main()
