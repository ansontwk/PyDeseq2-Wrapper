from utils import rna_analysis as analysis

#this version of pydeseq2 use OOP for a wrapper for pydeseq
def main():
    
    airway_set = analysis.rna_dataset("./data/airway_rawcounts.csv",
                                    "./data/airway_metadata.csv",
                                    ("control, treated"),
                                    treatmentheader = 'dex',
                                    name = 'airway')
    airway_set.initalize()
    airway_set.run_inference()
    airway_set.vst()
    airway_set.shrinklfc()
    
    genelist = ["ENSG00000152583", "ENSG00000179094", "ENSG00000116584", "ENSG00000189221", "ENSG00000120129", "ENSG00000148175"]
    
    analysis.plot_gene_compare(airway_set, 
                               gene = genelist, 
                               save_path = './plots/airway_compare_mult_t.png',
                               log_transform = True)
    analysis.plot_gene_compare(airway_set, 
                            gene = genelist, 
                            save_path = './plots/airway_compare_mult_f.png',
                            log_transform = False)
    analysis.plot_gene_compare(airway_set, 
                            gene = "ENSG00000152583",
                            save_path = './plots/airway_compare_single_t.png',
                            log_transform = True)    
    analysis.plot_gene_compare(airway_set, 
                            gene = "ENSG00000152583",
                            save_path = './plots/airway_compare_single_f.png',
                            log_transform = False)  
    
    analysis.plot_MA(airway_set,  
                     save_path="./plots/airway_MA_plot.png", 
                     s = 20)

    analysis.plot_PCA(airway_set,   
                      save_path='./plots/airway_PCA.png', 
                      s = 20)
    
    analysis.plot_volcano(airway_set, 
                          save_path='./plots/airway_volcano.png')
    
main()
