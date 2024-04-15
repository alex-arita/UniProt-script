# UniProt-script
This script was developed as part of the "Basic Concepts of Bioinformatics" course at the University of Valencia to fully automate the task of annotating data from databases with specific biological information of interest.

It takes individual files containing information about the experiment from which the data was obtained and the gene symbols for the genes to annotate. The structure of the files looks like this:

    #Alex Fernando Arita Arita
    #GSE212622
    #Homo sapiens | Mus musculus
    #non_infected vs infected
    Il21
    Tigit
    Il12rb1
    Cd6
    Cxcr3
    Ccl5
    Gbp2
    Icos
    ...

And returns new files in a "results" folder with the following type of information:
    #Alex Fernando Arita Arita
    #GSE212622
    #Homo sapiens | Mus musculus
    #non_infected vs infected
    "From"	"Entry"	"Entry Name"	"Reviewed"	"Protein names"	"Gene Names"	"Organism"	"Length"
    "Ackr4"	"Q9NPB9"	"ACKR4_HUMAN"	"reviewed"	"Atypical chemokine receptor 4 (C-C chemokine receptor type 11) (C-C CKR-11) (CC-CKR-11) (CCR-11) (CC chemokine receptor-like      1) (CCRL1) (CCX CKR)"	"ACKR4 CCBP2 CCR11 CCRL1 VSHK1"	"Homo sapiens (Human)"	350
    "Aif1"	"P55008"	"AIF1_HUMAN"	"reviewed"	"Allograft inflammatory factor 1 (AIF-1) (Ionized calcium-binding adapter molecule 1) (Protein G1)"	"AIF1 G1 IBA1"	"Homo          sapiens (Human)"	147
    "Apol6"	"Q9BWW8"	"APOL6_HUMAN"	"reviewed"	"Apolipoprotein L6 (Apolipoprotein L-VI) (ApoL-VI)"	"APOL6 UNQ3095/PRO21341"	"Homo sapiens (Human)"	343
    "Bcl3"	"P20749"	"BCL3_HUMAN"	"reviewed"	"B-cell lymphoma 3 protein (BCL-3) (Proto-oncogene BCL3)"	"BCL3 BCL4 D19S37"	"Homo sapiens (Human)"	454

The scripts is a modified version from Programmatic Access - ID Mapping by UniProt.
