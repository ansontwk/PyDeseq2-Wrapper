import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.decomposition import PCA
import seaborn as sns
plt.rcParams["font.family"] = "serif"

class rna_dataset():
    def __init__(self, 
                 raw_counts_path: str, 
                 metadata_path: str, 
                 columns: tuple, 
                 treatmentheader: str,
                 name: str | None = None, 
                 fix_raw: bool = True,
                 filter_counts: bool = True,
                 filter_threshold: int = 10):
        
        self.name = name #name of dataset
        self.raw_counts_path = raw_counts_path
        self.metadata_path = metadata_path
        
        #self.dataloader = rnaseq_dataloader(self.raw_counts_path, self.metadata_path, self.column)
        #self.raw_counts = raw_counts.apply(pd.to_numeric)
        #self.metadata = metadata
        self.column = str(columns)
        self.treatment_header = treatmentheader

        self._read_file()
        if fix_raw:
            self._fix_rawcounts()
        if filter_counts:
            self.filter_counts(threshold = filter_threshold)
            
        self.genes = list(self.raw_counts.columns)
        self.samples = list(self.raw_counts.index)
        
    def _read_file(self):
        #accepted format with 1 sample each row. Gene on columns
        self.raw_counts = pd.read_csv(self.raw_counts_path)
        self.metadata = pd.read_csv(self.metadata_path, index_col=0) 
        #ids must be row name   
        
    def filter_counts(self, threshold = 10):
        #filters gene counts with less than a threshold
        keep = self.raw_counts.sum(axis=0) >= threshold
        self.raw_counts = self.raw_counts.loc[:, keep]

    def initalize(self, ncpu = 10):
        #Setup analysis
        self.inference = DefaultInference(n_cpus=ncpu)
        design = f"~ {self.treatment_header}"
        self.dds = DeseqDataSet(counts = self.raw_counts,
                           metadata = self.metadata,
                           design = design,
                           refit_cooks= True,
                           inference= self.inference)
    
    def run_inference(self, set_norm: bool = True):
        self.dds.deseq2()
        header = [self.treatment_header]
        cols = self.column.split(",")
        for group in cols:
            header.append(group.strip())

        self.ds = DeseqStats(self.dds, contrast = header, inference = self.inference)
        self.ds.summary()
        
        if set_norm:
            self._norm_count()
        
    def _get_result(self):
        #deprecated
        print(self.ds.results_df)

    def vst(self, return_data: bool = False):
        self.dds.vst()
        vst = self.dds.layers["vst_counts"]
        vst_df = pd.DataFrame(vst, columns = self.genes, index = self.samples)
        self.vst_df = vst_df
        if return_data:
            return vst_df
    
    def _norm_count(self, return_data: bool = False):
        norm_counts = self.dds.layers['normed_counts']
        norm_counts_df = pd.DataFrame(norm_counts, columns = self.genes, index = self.samples)
        self.norm_count = norm_counts_df
        
        if return_data:
            return norm_counts_df
        
        #print(vst.shape)
        #normie = self.dds.layers["normed_counts"] #normalize
        #print(normie, normie.shape)
    
    def shrinklfc(self, col = 1):
        #for lfc shrinkage for post processing and visualization
        #call this to shrink coefficients
        self.ds.lfc_shrink(coeff=self.dds.obsm["design_matrix"].columns[col])
        #this directly overrides the original values in results df
        
    def check_shrink_lfc(self):
        #use this to check the design matrix 
        #in case you're not sure if the lfc shrinkage is correctly configured
        print(self.dds.obsm["design_matrix"].columns)
        print(self.dds.varm["LFC"])
    
    def _fix_rawcounts(self):
        #transposes the raw count matrix, and fix headers.
        self.raw_counts = remove_header(transposefit(self.raw_counts)).apply(pd.to_numeric)
    
'''class rnaseq_dataloader():
    #depracated. dataloader now incorporated into the RNA dataset class
    def __init__(self, raw_path, meta_path, condition_column):
        self.raw_path = str(raw_path)
        self.meta_path = str(meta_path) 
        self.condition_column = condition_column #accepts as list
    
    def read_file(self):
        #accepted format with 1 sample each row. Gene on columns
        self.raw_count = pd.read_csv(self.raw_path)
        self.metadata = pd.read_csv(self.meta_path, index_col=0) 
        #ids must be row name'''
    
def transposefit(pdframe):
    pdframe = pdframe.transpose()
    pdframe.columns = pdframe.iloc[0]
    pdframe = pdframe[1:]
    return pdframe

def remove_header(pdframe):
    #removes the top left corner item in columns
    pdframe = pdframe.rename_axis(None, axis = 1)
    return pdframe
        
def plot_MA(rna_dataset,
            colour: str = "#3C5488FF",
            logscale: bool = True, 
            save_path: str | None = None,
            dataset_name: str | None = None,
            **kwargs):
    #receives rna_dataset object to make a MA plot
    if rna_dataset.name is not None and dataset_name is None:
        dataset_name = rna_dataset.name
    colours = rna_dataset.ds.results_df['padj'].apply(lambda x: colour if x < 0.05 else "grey")
    
    fig, ax = plt.subplots(dpi=600)
    
    kwargs.setdefault('alpha', 0.5)
    kwargs.setdefault('s', 0.2)
    
    plt.scatter(x = rna_dataset.ds.results_df['baseMean'],
                y = rna_dataset.ds.results_df['log2FoldChange'],
                c = colours,
                **kwargs)
    
    ax.set_adjustable("datalim")
    
    if logscale:
        plt.xscale('log')
        
    plt.axhline(0, color="#DC0000FF", alpha=0.5, linestyle="--", zorder=3)    
    plt.xlabel('Mean of normalized counts')
    plt.ylabel('Log2 fold change')
    if dataset_name is not None:
        plt.title(f'M-A plot of {dataset_name}')
    else:
        plt.title('M-A plot')
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path, bbox_inches = 'tight')    
        plt.clf()
        
def plot_PCA(rna_dataset, 
             colour_palette = ('#3C5488FF', '#DC0000FF'),
             save_path: str | None = None,
             dataset_name: str | None = None,
             **kwargs):
    treatmentheader = rna_dataset.treatment_header
    if rna_dataset.name is not None and dataset_name is None:
        dataset_name = rna_dataset.name
    #plots PCA, only run this after running vst() function
    #accepts rna_dataset object
    
    kwargs.setdefault('alpha', 0.5)
    kwargs.setdefault('s', 0.2)
    
    vst = rna_dataset.dds.layers['vst_counts']
    X = np.asarray(vst)
    pca = PCA(n_components=2).fit_transform(X)
    conditions = rna_dataset.metadata[treatmentheader].values.tolist()
    condition_uniq = list(set(conditions))

    colours = [colour_palette[0] if i == conditions[0] else colour_palette[1] for i in conditions]
    fig, ax = plt.subplots(dpi=600)
    
    plt.scatter(pca[:,0], pca[:,1], c=colours, **kwargs)
    for i, sid in enumerate(rna_dataset.metadata.index):
        plt.text(pca[i,0], pca[i, 1], sid, fontsize = 8, wrap = True)
        
    ax.set_adjustable("datalim")
    plt.xlabel("PC1")
    plt.ylabel("PC2")

    handles = [Rectangle((0, 0), 1, 1, color=colour_palette[i]) for i in range(len(condition_uniq))]
    plt.legend(handles, 
               condition_uniq, 
               title='Condition', 
               bbox_to_anchor=(1.04, 0.5), 
               loc="center left", 
               borderaxespad=0)
    
    if dataset_name is not None:
        plt.title(f"VST-PCA by condition on {dataset_name}")
    else:
        plt.title("VST-PCA by condition")
    
    #plt.tight_layout()
    
    if save_path is not None:
        plt.savefig(save_path, bbox_inches = 'tight')
    else:
        plt.show()
    plt.clf()
    
def _maplog(x):
    return np.log(x + 0.5)

def _volcanolog(x):
    return -np.log(x)
 
def plot_gene_compare(rna_dataset,
                      gene: str | list,
                      normalized: bool = True,
                      log_transform: bool = False,
                      save_path: str | None = None,
                      **kwargs):
    treatmentheader = rna_dataset.treatment_header
    #plots gene comparison between two groups    
    counts = rna_dataset.norm_count[gene]
    conds = rna_dataset.metadata[treatmentheader]
    
    if log_transform:
        if type(gene) == str:
            counts = counts.map(_maplog)
        else:
            counts = counts.applymap(_maplog)
            
    df = pd.concat([counts, conds], axis = 1) 

    if type(gene) == str:
        #single gene
        sns.swarmplot(df, 
                    y = df[gene], 
                    x = df[treatmentheader],
                    orient = 'v',
                    color = '#3C5488FF')
        plt.xlabel('Condition')
        if log_transform:
            plt.ylabel('Log Normalized count')
        else:    
            plt.ylabel('Normalized Count')
        plt.title(f"{gene}")
    
    else:
        #multigene
        groups = list(set(df[treatmentheader].values.tolist()))
        df_melt = df.melt(id_vars = treatmentheader, var_name = 'Gene', value_name = 'Count', value_vars = gene)
        g = sns.FacetGrid(df_melt, col="Gene", col_wrap = 3, sharey = False)
        g.map(sns.swarmplot, treatmentheader, 'Count', color = '#3C5488FF', order = groups)
        
        if log_transform:
            g.set_axis_labels('Condition', 'Log Normalized count')
        else:
            g.set_axis_labels('Condition', 'Normalized count')
        
    if save_path is not None:
        plt.savefig(save_path, bbox_inches = 'tight', dpi = 600)
    else:
        plt.show()
    plt.clf()

def plot_volcano(rna_dataset,
                save_path: str | None = None,
                dataset_name: str | None = None,
                **kwargs):
    if rna_dataset.name is not None and dataset_name is None:
        dataset_name = rna_dataset.name
    kwargs.setdefault('alpha', 0.5)
    kwargs.setdefault('s', 0.45)
    fig, ax = plt.subplots(dpi=600)

    df = rna_dataset.ds.results_df[['log2FoldChange', 'pvalue', 'padj']]
    df['-log(pvalue)'] = df['pvalue'].map(_volcanolog)
    
    colours = [] # a bit hacky, but eh
    for _, row in df.iterrows():
        if row['log2FoldChange'] < 0 and row['padj'] <= 0.05:
            colours.append('#DC0000FF')
        elif row['log2FoldChange'] > 0 and row['padj'] <= 0.05:
            colours.append('#3C5488FF')
        else:
            colours.append('grey')    
            
    plt.scatter(x = df['log2FoldChange'], y = df['-log(pvalue)'], c = colours, **kwargs)
    if dataset_name is not None:
        plt.title(f"Volcano plot of {dataset_name}")
    else:
        plt.title(f"Volcano plot")
        
    plt.xlabel('logFoldChange')
    plt.ylabel('-log10(pvalue)')
    
    ax.set_adjustable("datalim")
    if save_path is not None:
        plt.savefig(save_path, bbox_inches = 'tight')
    else:
        plt.show()
    plt.clf