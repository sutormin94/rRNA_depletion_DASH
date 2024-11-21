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

# Chromosome name from which rRNA region is pulled.
Chromosome_name_spec="AE004091.2"

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

def read_average_write_wig(dict_of_files, coordinates_ar, chromosome_name_spec, output_path, dataset_name):
    
    # Read coverage depth data in wig format.
    dict_of_replicas={}
    
    for wig_name, wig_path in dict_of_files.items():
    
        print('Now is processing: ' + str(wig_path))
        wigin=open(wig_path, 'r')
        Dict_of_chromosomes_seqs={}
        NE_values=[]
        for line in wigin:
            line=line.rstrip().split(' ')
            if line[0] not in ['track', 'fixedStep']:
                NE_values.append(float(line[0]))
                
            else: 
                if 'chrom=' in line[1]:
                    if len(NE_values)>0:
                        Dict_of_chromosomes_seqs[Chrom_name]=np.array(NE_values)
                    Chrom_name=line[1].split('=')[1]
                    print(Chrom_name)
                    
                    NE_values=[]
                    
        Dict_of_chromosomes_seqs[Chrom_name]=np.array(NE_values)        
        wigin.close()
        
        dict_of_replicas[wig_name]=Dict_of_chromosomes_seqs
        
    # Normalise coverage depth data by total coverage.
    coef=10e6
    for wig_name, wig_data in dict_of_replicas.items():
        total_coverage=0
        for Chrom_name, Chrom_data in wig_data.items():
            total_coverage=np.sum(Chrom_data)
            dict_of_replicas[wig_name][Chrom_name]=(Chrom_data/total_coverage)*coef
        
    # Average normalized data.
    Aver_cov_depth_tracks_dict={}
    
    for Chrom_name in list(dict_of_replicas[list(dict_of_replicas.keys())[0]].keys()):
        
        Chrom_len=len(dict_of_replicas[list(dict_of_replicas.keys())[0]][Chrom_name])
        
        Aver_cov_depth_chrom_track=[]
        
        for i in range(Chrom_len):
            
            av_data_position=[]
            
            for wig_name, wig_data in dict_of_replicas.items():
            
                Chrom_position_data=wig_data[Chrom_name][i]
                av_data_position.append(Chrom_position_data)
                
            aver_cov_depth=np.mean(av_data_position)
            Aver_cov_depth_chrom_track.append(aver_cov_depth)
            
        Aver_cov_depth_tracks_dict[Chrom_name]=Aver_cov_depth_chrom_track
                     
    
    # Select region of interest.
    Aver_cov_depth_ROI=Aver_cov_depth_tracks_dict[chromosome_name_spec][coordinates_ar[0]:coordinates_ar[1]]

    #Write file with avaraged normalized data for specified ROI.
    average_out=open(os.path.join(output_path, f'{dataset_name}_rRNA.wig'), 'w')
    average_out.write(f'track type=wiggle_0 name="{dataset_name}" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom={chromosome_name_spec}_{dataset_name}_rRNA start=1 step=1\n')
    
    for i in range(len(Aver_cov_depth_ROI)):
        average_out.write(f'{Aver_cov_depth_ROI[i]}\n')
    average_out.close()
    
    return


#######
#Read reference genome and extract ROI.
#######

def read_extract_seq(ref_gen_path, coordinates_ar, chromosome_name_spec, dataset_name, output_path):
    
    outfile=open(os.path.join(output_path, f'{dataset_name}_rRNA.fasta'), 'w')
    
    for record in SeqIO.parse(ref_gen_path, 'fasta'):
        record_id=record.id
        record_seq=record.seq
        
        if chromosome_name_spec==record_id:
            
            record_id_ext=f'{chromosome_name_spec}_{dataset_name}_rRNA'
            record_seq_ROI=str(record_seq)[coordinates_ar[0]: coordinates_ar[1]]
            
            outfile.write(f'>{record_id_ext}\n{record_seq_ROI}')
            
            break
        
    outfile.close()
    
    return


def wrapper_function(dict_of_cov_depth, ref_gen_path, coordinates_ar, chromosome_name_spec, dataset_name, output_path):
    
    # Read, normalize, and average cov depth. Return norm av cov depth for ROI.
    read_average_write_wig(dict_of_cov_depth, coordinates_ar, chromosome_name_spec, output_path, dataset_name)
    
    # Read and extract sequence of ROI.
    read_extract_seq(ref_gen_path, coordinates_ar, chromosome_name_spec, dataset_name, output_path)
    
    return

wrapper_function(Dict_of_cov_depth, Ref_gen_path, Coordinates_ar, Chromosome_name_spec, Dataset_name, Output_path)
