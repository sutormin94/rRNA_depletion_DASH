###############################################
##Dmitry Sutormin, 2023##
##Select genomic region of interest##

# Select genomic region to which sgRNAs to be designed.
# Normalizes coverage depth by total coverage depth and averages normalized cov depth between replicates provided.
# Outputs normalized avreged cov depth for a region of interest (ROI).
###############################################

#######
# Packages to be imported.
#######

import os
import numpy as np
from Bio import SeqIO

#################
### Variables to be defined.
#################

# Dict of coverage depth profiles.
Dict_of_cov_depth={'1' : 'C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\PAO1\Genome_coverage_depth\P3_TdT.wig',
                   '2' : 'C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\PAO1\Genome_coverage_depth\P3_TS_TdT.wig'}

# Path to reference genome.
Ref_gen_path="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\PAO1\Reference_genome\Pseudomonas_aeruginosa_pao1_GCA_000006765.1_ASM676v1_genomic.fna"

# Coordinates of a region to be extracted. If targeting rRNA operon, choose one in + orientation.
Coordinates_ar=[722090, 727260]

# Output path.
Output_path="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\Guides_design_DASH\PAO1\\"

# Dataset name.
Dataset_name="P3"

Output_path=os.path.join(Output_path, 'Targeting_region')

if not os.path.isdir(Output_path):
    os.mkdir(Output_path)
    
    
#######
#Average wig files and write WIG.
#######

def read_average_write_wig(dict_of_files, coordinates_ar, output_path, dataset_name):
    
    # Read coverage depth data in wig format.
    dict_of_replicas={}
    
    for wig_name, wig_path in dict_of_files.items():
    
        print('Now is processing: ' + str(wig_path))
        wigin=open(wig_path, 'r')
        NE_values=[]
        for line in wigin:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                NE_values.append(float(line[0]))
        wigin.close()
        
        dict_of_replicas[wig_name]=np.array(NE_values)
        
    # Normalise coverage depth data by total coverage.
    coef=10e6
    for wig_name, wig_data in dict_of_replicas.items():
        total_coverage=np.sum(wig_data)
        dict_of_replicas[wig_name]=(wig_data/total_coverage)*coef
        
    # Average normalized data.
    Aver_cov_depth_track=[]
    
    for i in range(len(dict_of_replicas[list(dict_of_replicas.keys())[0]])):
        av_data_position=[]
        for replica_name, replica_data in dict_of_replicas.items():
            av_data_position.append(replica_data[i])
        aver_cov_depth=np.mean(av_data_position)
        
        Aver_cov_depth_track.append(aver_cov_depth)
    
    # Select region of interest.
    Aver_cov_depth_ROI=Aver_cov_depth_track[coordinates_ar[0]:coordinates_ar[1]]

    #Write file with avaraged normalized data for specified ROI.
    average_out=open(os.path.join(output_path, f'{dataset_name}_rRNA.wig'), 'w')
    average_out.write(f'track type=wiggle_0 name="{dataset_name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom={dataset_name}_rRNA start=1 step=1\n')
    
    for i in range(len(Aver_cov_depth_ROI)):
        average_out.write(f'{Aver_cov_depth_ROI[i]}\n')
    average_out.close()
    
    return


#######
#Read reference genome and extract ROI.
#######

def read_extract_seq(ref_gen_path, coordinates_ar, dataset_name, output_path):
    
    for record in SeqIO.parse(ref_gen_path, 'fasta'):
        record_id=f'{dataset_name}_rRNA'
        record_seq=str(record.seq)[coordinates_ar[0]: coordinates_ar[1]]
        
    outfile=open(os.path.join(output_path, f'{dataset_name}_rRNA.fasta'), 'w')
    
    outfile.write(f'>{record_id}\n{record_seq}')
    
    outfile.close()
    
    return

    
    
def wrapper_function(dict_of_cov_depth, ref_gen_path, coordinates_ar, dataset_name, output_path):
    
    # Read, normalize, and average cov depth. Return norm av cov depth for ROI.
    read_average_write_wig(dict_of_cov_depth, coordinates_ar, output_path, dataset_name)
    
    # Read and extract sequence of ROI.
    red_extract_seq(ref_gen_path, coordinates_ar, dataset_name, output_path)
    
    return

wrapper_function(Dict_of_cov_depth, Ref_gen_path, Coordinates_ar, Dataset_name, Output_path)
