# MINDPROT
MINDPROT: Markov Inside for Drugs and Proteins

Background:

MINDPRROT application based on the MARCH-INSIDE: Markov Chain Invariants for Networks Simulation and Design MARCH-INSIDE (MI) https://github.com/glezdiazh/MARCH-INSIDE. MI is a well-known method introduced by Prof. Humbert G Díaz (Gonzaléz-Díaz et al.) as early as 2002 for the calculation of Markov Invariants (Moments, Shanon entropies, Mean Markov values) of molecular graphs and complext netxorks using a Markov chain stchastic approach. In case you want to develop new collaborations, applications, etc. related to MI algorith please do not hesitate to contact us at: https://www.linkedin.com/in/humbertgdiaz/

MINDPROT Specificities:

MInD-Prot is a Python/wxPython application for the calculation of the Mean properties Markov indices for drugs and proteins.
The application can calculate averaged indices for proteins and drugs by using input classes (for drugs and proteins) or PDB origin (only for proteins). In addition, MInD-Prot can generate mixed indices for protein, drug or protein-drug pairs. If it is necessary, the tool can generate random negative pairs for protein and protein-drug pairs. Additional information can be obtained such as the PDB headers for proteins. These numbers that characterise each protein/drug or pair protein-drug is used for building QSAR/QPDRclassification models. 

MINDPROT Authors / Contributions:

Prof. Cristian Robert Munteanu (MINDPROT Software rograming, AI/ML applications, Co-author of papers, https://github.com/muntisa)

Prof. Humbert G. Díaz (MINDPROT algorithm and software design, AI/ML applications, Co-author of papers, https://github.com/glezdiazh),


How to use MINDPROT:
In the main windows you can choose the calculation parameters and the input/output parameters.
The main GUI is divided in the following parts:
PROTEINS
Files parameters
- input PDBchain list file as PDB/PDBchain
- protein simple result file
- the PDB folder local database; if the PDB is missing, it will be downloaded from the online PDB databank; if the PDB don't exists, an error messages will be printed

Alpha Carbon Network parameters
- Cutoff, Roff, Ron parameters for deciding if two amino acid carbon alpha are linked
- Protein Orbital Region Limits in % as core, inner, middle, outer
- By Chain calculation: if you enable this control, the tool will consider the entire protein and will calculate all the chains even if you have PDBchains in the input; if this is disabled, the calculation will consider each PDBchain in the list (if the chain info is missing, the entire protein network will be constructed)

Averaged Indices (PROT_ClassAvgs.txt)
- by PDB header information such as head, expression_system, expression_system_taxid, name, chain, organism_scientific, molecule, expression_system_vector_type, ec, organism_common, expression_system_plasmid, engineered, expression_system_strain, cell_line, cellular_location, gene, organism_taxid
- by input classe in the PDBchain list file as PDBchain[tab]Class

Full header information output (PROT_FullHeaderRes.txt)
- get the full info from PDBs and add them to the simplest result (one column for each header field)

Protein PAIRS
- using the similarity PDB of the chains as positive pairs and generating the negative ones until X times the positive pairs
- as alternative: use the activity file (default: ProtPairActivity.txt) with the PDB1[tab]PDB2[tab]Class

DRUGS
Files parameters
- SMILE list file (Drug Name[tab]SMILE formula)
- drug simple result file
- Averaged results by input classes (DRUG_ClassAvg.txt)

Averaged result (DRUG_ClassAvgs.txt)
- using the input classes from the SMILE list file (Drug Name[tab]SMILE formula[tab]Class)

Drug PAIRS
- using always the input activity file as DrugName1[tab]DrugName2[tab]Class; the identificacion of the drugs in the results is made by the drug name, not by the SMILE formula


PROTEIN-DRUG PAIRS
- using always an input activity file as PDBChain[tab]DrugName[tab]Activity
- if there is only one type of class (positive cases), it can be generated random protein-drug pairs until X times the positives
- you can calculate this type of pairs only if both PROTEIN and DRUG calculations are enabled

Notes:
You can create/edit and browse all the files directly from the interface by using the native NotePad from Windows.
All the options have default values in order to perform minimum of calculations. The processing of the network can be seen in a console window. If you will close it, all the application windows will close too.

If MI-Prot will raise any error, you can see it in the same console and please send us in order to fix it. Markov Indices from MInD-Prot Before to calculate the indices, the connectivity matrix is Markov normalized (the matrix element is divided by the maximum in its matrix row) resulting a matrix of node probabilities (P). In the second step, these P matrix is raised to the power k=5, resulting k matrices (Pk), the input for all the indice calculations. The node Markov centralities are calculated as the difference between the indices of the original matrix and the indices without the correspondent node. The MCs are based on the Mean Properties (MP) addapted Markov indices.
