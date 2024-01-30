###############################################
##Dmitry Sutormin, 2023##
##Rarefy sgRNAs##

# Takes sgRNAs designed by design_grnas.py (from the DASH paper).
# Takes coverage depth of a region targeted by scRNAs.
# Rarefies sgRNAs to a desired number by treating coverage depth as a probability density function to keep a sgRNA.
# Outputs rarefied guides as a csv table and guides coordinates in a bed file.
###############################################

#######
# Packages to be imported.
#######

import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

#################
### Variables to be defined.
#################

# Path to coverage depth profile (prepared by the Select_genomic_ragion.py script).
Cov_depth_profile="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\Guides_design_DASH\E_coli\Targeting_region\C2_subset_rRNA.wig"

# Path to a sequence of ROI (prepared by the Select_genomic_ragion.py script).
Ref_seq_path="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\Guides_design_DASH\E_coli\Targeting_region\C2_subset_rRNA.fasta"

# Path to guides compendium (prepared by the design_grnas.py script).
Guides_compend_path="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\Guides_design_DASH\E_coli\oligos_Quoc.csv"

# Handle sequence (exactly as was provided in design_grnas.py script).
Handle_seq="GTTTTAGAGCTAGA"

# Output_path.
Output_path="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\Guides_design_DASH\E_coli\\"

# Desired number of guides.
Total_num_guides=300

# Dataset name.
Dataset_name="E_coli_weighted"

def read_wig(wig_path):
    
    wigin=open(wig_path, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
        
    return np.array(NE_values)


def read_ref_seq(ref_seq_path):
    
    for record in SeqIO.parse(ref_seq_path, 'fasta'):
        record_id=record.id
        record_seq=str(record.seq)
        
    return record_id, record_seq


def read_guides(guides_compend_path, handle_seq):
    
    Guides_dict={}
    
    filein=open(guides_compend_path, 'r')
    
    for line in filein:
        line=line.rstrip().split(',')
        guide_id=line[0]
        guide_full_seq=line[1]
        spacer_seq=line[1][-len(handle_seq)-20:-len(handle_seq)]
        Guides_dict[guide_id]=[spacer_seq, guide_full_seq]
    
    filein.close()
    
    print(f'Initial number of guides: {len(Guides_dict)}')
    
    return Guides_dict


def rarefy_guides(cov_depth_track, ref_seq_seq, Guides_dict, total_num_guides):
    
    # Identify coordinates of guides.
    Coords_dict={}
    lost_counter=0
    
    for guide_id, guide_data in Guides_dict.items():
        
        spacer_seq=guide_data[0]
        
        if spacer_seq in ref_seq_seq:
            start=ref_seq_seq.find(spacer_seq)
            end=start+len(spacer_seq)
            Coords_dict[guide_id]=[start, end]
        else:
            spacer_seq_Seq=Seq(spacer_seq)
            spacer_seq_rc=str(spacer_seq_Seq.reverse_complement())
            if spacer_seq_rc in ref_seq_seq:
                start=ref_seq_seq.find(spacer_seq_rc)
                end=start+len(spacer_seq_rc)   
                Coords_dict[guide_id]=[start, end]
            else:
                lost_counter+=1
        
    print(f'Number of lost guides: {lost_counter}')
    print(f'Number of validated guides: {len(Coords_dict)}')
    
    # Get coverage depth for spacers.
    guides_cov_depth_ar=[]
    Cov_dict={}
    
    for guide_id, guide_coords in Coords_dict.items():
        
        guide_cov_depth=np.mean(cov_depth_track[guide_coords[0]:guide_coords[1]])
        guides_cov_depth_ar.append(guide_cov_depth)
        Cov_dict[guide_id]=guide_cov_depth
        
    Mean_guide_cov_depth=np.mean(guides_cov_depth_ar)
    print(f'Mean coverage depth of guides: {Mean_guide_cov_depth}')
    
    # Rarefy guides.
    Final_guides_dict={}
    for guide_id in Coords_dict.keys():
        
        guide_cov_depth=Cov_dict[guide_id]
        guide_full_seq=Guides_dict[guide_id][1]
        guide_prob=(float(total_num_guides)/len(Coords_dict))*(guide_cov_depth/Mean_guide_cov_depth)
        
        if guide_prob>1:
            Final_guides_dict[guide_id]=[guide_full_seq, Coords_dict[guide_id][0], Coords_dict[guide_id][1]]
        else:
            rand_value=float(np.random.rand(1,1))
            if guide_prob>rand_value:
                Final_guides_dict[guide_id]=[guide_full_seq, Coords_dict[guide_id][0], Coords_dict[guide_id][1]]
    
    print(f'Number of rarefied guides: {len(Final_guides_dict)}')
    
    return Final_guides_dict


def write_final_guides(dataset_name, Final_guides_dict, total_num_guides, ref_seq_id, output_path):
    
    # Write CSV file with selected guides.
    guides_csv=open(os.path.join(output_path, f'{dataset_name}_rarefied_{total_num_guides}_{len(Final_guides_dict)}.csv'), 'w')
    
    guides_csv.write(f'Guide ID\tPool name\tSequence\n')
    
    for guide_id, guide_info in Final_guides_dict.items():
        
        pool_name=guide_id.split('_')[0]
        guide_full_seq=guide_info[0]
        
        guides_csv.write(f'{guide_id}\t{pool_name}\t{guide_full_seq}\n')
    
    guides_csv.close()
    
    # Write BED file with coordinates of selected guides.
    guides_bed=open(os.path.join(output_path, f'{dataset_name}_rarefied_{total_num_guides}_{len(Final_guides_dict)}.bed'), 'w')
    
    for guide_id, guide_info in Final_guides_dict.items():
        
        guide_coord_start=guide_info[1]
        guide_coord_end=guide_info[2]
        
        guides_bed.write(f'{ref_seq_id}\t{guide_coord_start}\t{guide_coord_end}\t{guide_id}\n')    
    
    guides_bed.close()
    
    return


def wrapper_function(dataset_name, cov_depth_profile, ref_seq_path, guides_compend_path, handle_seq, total_num_guides, output_path):
    
    # Read coverage depth profile.
    cov_depth_track=read_wig(cov_depth_profile)
    
    # Read referense sequence.
    ref_seq_id, ref_seq_seq=read_ref_seq(ref_seq_path)
    
    # Read initial guides compendium to be rarefied.
    Guides_dict=read_guides(guides_compend_path, handle_seq)
    
    # Rarefy guides.
    Final_guides_dict=rarefy_guides(cov_depth_track, ref_seq_seq, Guides_dict, total_num_guides)
    
    # Write final data: guides as a csv table; guides coordinates as a bed file.
    write_final_guides(dataset_name, Final_guides_dict, total_num_guides, ref_seq_id, output_path)
    
    return

wrapper_function(Dataset_name, Cov_depth_profile, Ref_seq_path, Guides_compend_path, Handle_seq, Total_num_guides, Output_path)