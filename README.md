## ROHs  
### scripts to generate binary ROH files for PLINK .hom out  
### scripts to plot ROHs per sample (animal)


## A) Generate ROH binary file for PLINK output .hom 
### Using the script "ROH_binary.R"  

#### Input files (5 argument used)  
    - rohsummary ==== name of PLINK '.hom' file  
    - mapfile    ==== mapfile used to generate the ROH files  
    - animlist   ==== file containing all animals used in generating ROH files  
    - chrNumber  ==== just a single chromosome number is allowed  
    - outputname ==== outputname

#### Output files (2 files are produced)
    - "*_binary.roh"  ==== This contains the ROH binary file  
    - ".snplist"      ==== SNP list  


#### see 'example_roh_binary.R' for further details




### B) Plot ROH per sample using the binary ROH file generated from the "ROH_binary.R" script  
### Using the script "plot_ROHs.R"  

#### Input files (5 argument used)  
    - rohbinaryfile  ==== name of binary roh file  
    - mapfile    ==== mapfile used to generate the ROH files  
    - chrNumber  ==== just a single chromosome number is allowed  
    - anim   ==== the row numbers of the animals you want to plot  
    - colourRohs ==== do you want to colour the ROHs or the non-ROHs if YES, the rohs are coloured, if NO then non-ROHs are coloured  
    - title ==== title of the graph  
    - outname ==== output name of the graph  

#### Output files (1 files are produced)  
    - The graph with the outname - ".tiff"  

#### see 'example_plot_ROH.R' for further details  


