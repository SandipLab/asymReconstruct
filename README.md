# asymReconstruct
The asymReconstruct is a workflow first introduced in [Asymmetric gating of a homopentameric ion channel GLIC revealed by cryo-EM](https://www.pnas.org/doi/10.1073/pnas.2512811122#executive-summary-abstract), reported by Li et al., where it was successfully used to decipher the asymmetric conformational changes during the activation of a homopentameric proton-gated ion channel GLIC.    
  
The workflow is primarily based on cryo-EM data processing of single particles with cyclic pseudosymmetry, using the [RELION](https://github.com/3dem/relion) software,  with two stand-alone shell scripts to perform key steps, namely ```configurations_stats.sh``` and ```align.sh```.

A general demonstration of asymReconstruct workflow is described below using publicly available data.

## Case study: Reconstruction of asymmetric conformations of a homopentameric ligand-gated ion channel 5-HT3AR
Processing the restacked particles from EMPIAR-10381, including symmetry expansion, focused 3D classification, particle reorientation, and reconstruction/local refinement.
### Introduction
Symmetry and pseudosymmetry are common among various biomacromolecules. One of the prevalent types of symmetry is cyclic symmetry, denoted as Cn, where n is the number of subunits. Innumerable widely studied ion channels exhibit such symmetry, providing significant advantages for the cryo-EM reconstructions. However, conformational changes of homomeric complexes do not necessarily occur synchronously among all subunits, potentially creating asymmetric intermediate conformations. When such asymmetry is prominent, the cryo-EM reconstruction can be achieved simply without symmetry indicator, such as shown by [Matthies et al.](https://pubmed.ncbi.nlm.nih.gov/26871634/). However, comprehensive analysis of asymmetric conformational changes remains difficult due to the challenges in precise reconstruction of subtly asymmetric conformations present in a heterogeneous cryo-EM dataset.

To identify asymmetric conformations, 3D classifications are usually applied to symmetry-expanded particles using a single-subunit mask, such as demonstrated by [Roh et al.](https://pubmed.ncbi.nlm.nih.gov/28710336/), [Zhang et al.](https://pubmed.ncbi.nlm.nih.gov/33594077/), [Lee et al.](https://pubmed.ncbi.nlm.nih.gov/36805660/), and many others. However, the reconstruction of one or several asymmetric conformations remains challenging. 

Recently, we reported the asymmetric activation of a bacterial homopentameric proton-gated ion channel, GLIC, showing that all five subunits undergo asynchronous rotations and trigger expansion of a characteristic loop during gating (PMID: [41129221](https://pubmed.ncbi.nlm.nih.gov/41129221/)). On top of subunit-focused 3D classifications, we have used a simple strategy to reorient asymmetric particles which facilitate precise alignment for high-resolution 3D reconstructions.

While the details of asymmetric reconstructions of GLIC have been described in the publication and raw data are deposited at EMPIAR-13083 (TBR), we further demonstrate the general applicability of our method in the reconstructions of asymmetric conformations from pseudosymmetric single-particle datasets with cyclic symmetry, using the deposited dataset of a homopentameric serotonin receptor (EMPIAR-10381). For this dataset, the authors have reported two conformations of subunits according to focused 3D classifications, namely “more” (M, blue) and “less” (L, orange) open conformations (figure below). Backtracking these subunits to original particles identified all eight possible pentameric combinations or configurations (figure below). The cryo-EM structure of LMLML conformation was reported by the authors.
<p align="center">
  <img width="655" height="373" alt="fig1" src="https://github.com/user-attachments/assets/91c4a916-0d4e-46c9-9fc8-149e0ea2cb62" />
</p>  
  
### Requirements
RELION version 3 or newer. 

### 1. Particle stack import and reconstruction  
The files for the stacked final particles corresponding to serotonin-bound saposin-nanodisc-embedded 5-HT3AR can be downloaded from EMPIAR-[10381](https://www.ebi.ac.uk/empiar/EMPIAR-10381/), namely “
$${\color{red}251k\\_particles.mrcs}$$
” and “
$${\color{red}251k\\_particles.star}$$
” (figure below).
<p align="center">
  <img width="285" height="156" alt="fig2_screenshot_trim" src="https://github.com/user-attachments/assets/15cc480c-a7f2-4c6c-bf4f-96393a52515f" />
</p>

  
Notice that the image names inside the star file refer to polish output (highlighted in figure below).   
  
<img width="1503" height="108" alt="fig3" src="https://github.com/user-attachments/assets/c97bd9eb-c668-454e-8ced-f53f45629844" />      
  
We renamed that column inside the star file according to 
$${\color{red}251k\\_particles.mrcs}$$
using command below. The output file name is called as “
$${\color{red}251k\\_particles.cor.star}$$
”. Alternatively, you can rename that column using any advanced text editor (such as Geany).   

```awk '{if(NF<28)print;if(NF==28)print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, NR-33  "@251k_particles.mrcs", $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28}' 251k_particles.star > 251k_particles.cor.star```  
  
After renaming the star file column looks as below.  

<img width="1501" height="87" alt="fig4" src="https://github.com/user-attachments/assets/257f3a99-9afb-4326-b916-837bf573c1b0" />  

  
In RELION, reconstruct a 3D map using 
$${\color{red}251k\\_particles.cor.star}$$
using the alignment information already present in the metadata (this step along with previous steps can be ignored if you are processing your own data from scratch). The purpose of this step is to generate 
$${\color{red}initial\\ map}$$
for subsequent steps.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Input images STAR file:** 
$${\color{red}251k\\_particles.cor.star}$$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Reference map:** *emd-10692.mrc*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Symmetry:** *C5 (pentameric channel)*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Use GPU accelerations?:** *No*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Additional arguments:** *--skip_align*  

  
Alternatively, this part can be done in [cryoSPARC](https://cryosparc.com/). 
### 2. 3D classification focused on single subunit in RELION
To perform focused 3D classification in RELION, the following files are required:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
$${\color{red}Initial\\ map}$$
(from the last section)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
$${\color{blue}A\\ single\\ subunit\\ tight\\ mask}$$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
$${\color{blue}A\\ single\\ subunit\\ initial\\ volume}$$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
$${\color{blue}A\\ star\\ file\\ corresponding\\ to\\ symmetry\\ expanded\\ particles}$$  

To generate a single-subunit tight mask, different methods can be used. The standard approach is to gradually erase the densities of four subunits in ChimeraX and use the resulted single-subunit density to generate a mask in either Relion or CryoSPARC, with appropriate parameters including viewing threshold, dilation radius and soft-padding radius. All these parameters should be tested and optimized to ensure effective focused classification.

However, for this documented dataset, an atomic model is readily available in Protein Data Bank, with PDB ID 6Y5A. In such case, the single subunit mask can be easily prepared from the atomic model, by first overlaying the model to the 
$${\color{red}initial\\ map}$$ reconstructed as in the last step (figure below, left), and then generate a single-subunit 3D volume from the pdb using the following command in ChimeraX:

```
molmap #2/A 3.5 onGrid #1
```

where the “#2/A” denotes chain A of model 6Y5A, “3.5” indicates desired resolution to generate the volume, the “onGrid #1” ensures the box size of generated volume is consistent with that of the full-protein initial map, which is denoted by “#1”. The resulting generated volume is shown in the figure below (middle). #1 and #2 denote the map or model ID in ChimeraX.

Subsequently, the 
$${\color{blue}single\\ subunit\\ tight\\ mask}$$ can be generated in either Relion or cryoSPARC. A representative single-subunit mask is shown below in grey volume (right). 

<p align="center">
  <img width="631" height="289" alt="fig5" src="https://github.com/user-attachments/assets/c25106dd-8511-43c7-bced-7ab581313a9b" />
</p>

The 
$${\color{blue}single\\ subunit\\ initial\\ volume}$$ can be generated by multiplying the mask with the 3D reconstruction using the following command in ChimeraX:
```
vop multiply #1 #2 onGrid #1
```
where the “#1” points to the 3D reconstruction and “#2” points to the single-subunit mask, the “onGrid” option ensures that the box size of the volume product is consistent with the 3D reconstruction. The 
$${\color{blue}symmetry\\ expanded\\ star\\ file}$$ will be generated using the following command:
```
relion_particle_symmetry_expand --i 251k_particles.cor.star --o 251k_particles.cor.expanded.star --sym C5
```

With all input files prepared, the focused 3D classification can be performed while specifying the following parameters:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Input images STAR file:** 
$${\color{blue}251k\\_particles.cor.expanded.star}$$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Reference map:** <name of 
$${\color{blue}single\\ subunit\\ initial\\ volume}$$
\>  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Reference mask:** <name of 
$${\color{blue}single\\ subunit\\ tight\\ mask}$$
\>   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Symmetry:** C1  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Number of classes:** 2  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Regularization parameter Τ:** 16  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Perform image alignment?:** No  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Use GPU acceleration?:** No  

Importantly, different values and combinations of class numbers and regularization parameters should be tested for a new project to optimize the quality of 3D classifications and identify all possible conformations. For this demonstration, a class number of “2” is used based on the existing knowledge of M and L conformations shown by the authors of this dataset. 

Through this focused 3D classification, all subunits from 251k particles have been split into two classes corresponding to M (class 1) and L (class 2), consisting of 36.6% and 63.4% of the population, respectively (figure below). The major differences between these two conformations are seen at the transmembrane helices and the intracellular MX helix. The run_it025.star file (suppose 25 iterations) from the 3D classification is used for subsequent analysis.

<p align="center">
  <img width="790" height="345" alt="fig6" src="https://github.com/user-attachments/assets/d1f9e820-2a1b-47df-8e79-55978bb67d6c" />
</p>

### 3. Configurations statistics using ```configurations_stats.sh```
As has been demonstrated by multiple studies, it is feasible to perform particle configuration analysis based on the classification result. Here, we provide a simple and generalized script for configurations statistics to assist this workflow, namely configurations_stats.sh.

An output star file from the focused 3D classification job is required as the input, which contains subunit conformation assignment signified by the assigned class number. Thus, the statistics can be run with the following command:

```./configurations_stats.sh <name_ run_it025.star> ML 5```

where M indicates class 1 and L indicates class 2. “5” specifies the number of subunits for each particle, in this case, referring to a pentamer. More information is available in the documentation within the script. The output shows all existing configurations with subunits arranged in a clockwise order viewing against Z axis (figure below). Notably, the five letters in each configuration name constitute a cyclic list, thus, “MLLLL” and “LLLLM” are rotationally equivalent and refer to the same configuration.

<p align="center">
  <img width="651" height="108" alt="fig7" src="https://github.com/user-attachments/assets/3c79f369-10e3-48e9-ad56-024645de2837" />
</p>

### 4. Particle reorientation using ```align.sh```
Following the identification of different asymmetric configurations, subsets of particles belonging to each configuration can be isolated into separate star files and thus discussed separately. Taking MLLLL as an example. Isolating the particles with four subunits assigned to class 2 (L) and one subunit assigned to class 1 (M) during focused 3D classification results in a conformationally homogeneous dataset. However, conventional 3D refinements with symmetry-relaxation do not necessarily achieve optimized alignment among subunits, depending on the extent of asymmetric features. Consequently, the M subunit in each particle image may not align to the same subunit in the resulting 3D volume.

A convenient solution is to reorient particles to ensure subunit correspondence based on the run_it025.star file. The align.sh script performs such angular rearrangement by using the following command:

```./align.sh <name_run_it025.star> ML 5 MLLLL <name_output.star>```

where M and L indicate class 1 and class 2, “5” indicates pentamer, “MLLLL” specifies the configuration to reorient, and output provided by the user. An output star file name is needed. More information is available in the documentation within the script.

### 5. Reconstruction or local refinement
The output star file from align.sh is readily usable for either reconstruction or refinements with local angular search. 

In RELION, use the low-pass-filtered C5 map as the 
$${\color{red}initial\\ map}$$ and impose C1 symmetry. For simple reconstruction, specify:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Use GPU acceleration?:** *No*  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Additional arguments:** *--skip_align*  

For refinement with local angular search, specify:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Additional arguments:** *--sigma_ang 4*  

Note that the value for angular search range should be tested for a new project to ensure that there is no significant divergence of rotational angles after the refinement. For a pentameric protein, the change in *_rlnAngleRot* for each particle should be less than ±72 degrees.

To speed up the process, [cryoSPARC](https://cryosparc.com/) can be used. In cryoSPARC, simple reconstruction can be achieved by “Reconstruct Only” job as described above. And a refinement with local angular search can be performed further using the “Local Refinement” job. Users are advised to randomly split data into even halves for validation and resolution estimation.

Using cryoSPARC’s reconstruction, we have reconstructed seven conformations from the original 251k particles, with resolutions ranging from 2.96 to 3.66 Å (figure below).

<p align="center">
  <img width="752" height="388" alt="fig8" src="https://github.com/user-attachments/assets/34ee45d0-087a-428d-a9d9-8be5c14dbda1" />
</p>

### 6. Asymmetric features observed in the reconstructed maps
Further observations of the reconstructed maps proved the precise alignment of asymmetric features, displayed here as the extent of MX helix expansion (figure below).

<p align="center">
  <img width="752" height="454" alt="fig9" src="https://github.com/user-attachments/assets/f5a7fcdc-d0b8-4d74-988d-2bcd9b51a449" />
</p>










