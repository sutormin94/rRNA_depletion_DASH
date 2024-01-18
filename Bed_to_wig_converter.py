###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#Converts bed-like file(s) to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the working directory.
PWD="C:\\Users\dsutormi\D_Sutormin\Science\\rRNA_depletion_project\PAO1\Genome_coverage_depth\\"
#Path to the input file
filein_path_dict={'1' :  PWD + "P2_subsample_cov_depth.bed",  
                  '2' :  PWD + "P3_TdT_CKDL230034333-1A_HGNKFDSX7_L3_sorted.bed", 
                  '3' :  PWD + "P3_TS_TdT_CKDL230034333-1A_HGNKFDSX7_L3_sorted.bed", 
                  }

#Path to the output file.
fileout_path_dict={'1' :  PWD + "P2_subsample_cov_depth.wig", 
                   '2' :  PWD + "P3_TdT.wig", 
                   '3' :  PWD + "P3_TS_TdT.wig", 
                    }

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' :  "P2_subsample", 
           '2' :  "P3_TdT",
           '3' :  "P3_TS_TdT",
           }      

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)