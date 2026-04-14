# -*- coding=utf-8 -*-
# 2022.0114
# @zlk
import sys
import multiprocessing
import argparse
import numpy as np
import gzip

def read_vcf_data(infile,p1,p2,b1,b2):
    def get_geno(geno_str):
        # 修复版：增加容错处理
        try:
            tmp = geno_str.split(':')[1]
            if ',' not in tmp: return '0|0'
            parts = tmp.split(',')
            return parts[0] + '|' + parts[1]
        except:
            return '0|0'

    out_dict={}
    # 修复版：自动识别 gzip
    if infile.endswith('.gz'):
        vcf_file = gzip.open(infile, 'rt')
    else:
        vcf_file = open(infile, 'r')
        
    out_table_file=open(infile+'.table','w')
    for line in vcf_file:
        if line.startswith('#'):
            continue
        line_list=line.strip().split()
        # 增加长度检查防止越界
        if len(line_list) < max(p1, p2, b1, b2) + 1:
            continue
            
        if ',' in line_list[3] or ',' in line_list[4]:
            continue
        # filter genotrype by parents
        # for F2
        if line_list[p1].split(':')[0]==line_list[p2].split(':')[0] or line_list[p2].startswith('./.') or line_list[p1].startswith('./.') :
            continue
        # for F2
        elif line_list[b2].startswith('./.') or line_list[b1].startswith('./.'):
            continue
        elif line_list[p2].startswith('0/1') or line_list[p2].startswith('0|1') or line_list[p1].startswith('0/1') or line_list[p1].startswith('0|1'):
            continue
        p1_geno=get_geno(line_list[p1])
        if line_list[p2].startswith('./.') :
          p2_geno='0|0'
        else:
          p2_geno=get_geno(line_list[p2])
        # p2_geno=get_geno(line_list[p2])
        b1_geno=get_geno(line_list[b1])
        b2_geno=get_geno(line_list[b2])
        out_table_file.write(line_list[0]+'\t'+line_list[1]+'\t'+line_list[3]+'\t'+line_list[4]+'\t'+p1_geno+'\t'+p2_geno+'\t'+b1_geno+'\t'+b2_geno+'\n')
        if line_list[0] not in out_dict.keys():
            out_dict[line_list[0]]=[]
        out_dict[line_list[0]].append([int(line_list[1]),line_list[3],line_list[4],p1_geno,p2_geno,b1_geno,b2_geno])
    vcf_file.close()
    return out_dict

def read_table_file(infile):
    table_file=open(infile,'r')
    out_dict={}
    for line in table_file:
        line_list=line.strip().split('\t')
        if line_list[0] not in out_dict.keys():
            out_dict[line_list[0]]=[]
        ad_list=[int(line_list[1])]+line_list[2:]
        
        out_dict[line_list[0]].append(ad_list)
    return out_dict

def pre_binom_test(rep):
    out_dict={}
    for n in range(3,4000,1):
        num=np.percentile(np.random.binomial(n, 0.5, size=rep),99.99)
        out_dict[n]=num
    return out_dict

def filter_data(key,chr_list,dep_lim,dep_high,dep_lim_pool,dep_high_pool,binom_dict,snp_index,ED):
    out_snp_list=[]
    def get_dep(geno):
        # print(geno)
        try:
            dep=int(geno.split('|')[0])+int(geno.split('|')[1])
            geno=[int(geno.split('|')[0]),int(geno.split('|')[1])]
            return dep,geno
        except:
            return 0, [0,0]

    def binom_index(b1_dep,b2_dep,binom_dict):
        num1=binom_dict.get(int(b1_dep), 0)
        num2=binom_dict.get(int(b2_dep), 0)
        if b1_dep == 0 or b2_dep == 0: return 0
        
        index1=(num1/b1_dep-(b2_dep-num2)/b2_dep)**2
        
        index2=((b1_dep-num1)/b1_dep-num2/b2_dep)**2
        # print(num1,num2,b1_dep,b2_dep)
        # print(index1,index2)
        return index1+index2
    def binom_index2(b1_dep,b2_dep,binom_dict):
        num1=binom_dict.get(int(b1_dep), 0)
        num2=binom_dict.get(int(b2_dep), 0)
        if b1_dep == 0 or b2_dep == 0: return 0
        index1=(num1/b1_dep)
        index2=((b2_dep-num2)/b2_dep)
        return abs(index1-index2)
    for vcf_list in chr_list:
        p1_dep,p1_geno=get_dep(vcf_list[3])
        p2_dep,p2_geno=get_dep(vcf_list[4])
        b1_dep,b1_geno=get_dep(vcf_list[5])
        b2_dep,b2_geno=get_dep(vcf_list[6])
        # filter by depth
        if p1_dep<dep_lim or p1_dep>dep_high or p2_dep<dep_lim or p2_dep>dep_high or b1_dep<dep_lim_pool or b1_dep>dep_high_pool or b2_dep<dep_lim_pool or b2_dep>dep_high_pool:
            continue

        
        if p2_geno[1]>p2_geno[0]:
            p1_num,p2_num=0,1
        else:
            p1_num,p2_num=1,0
        if len(vcf_list[1])>1 or  len(vcf_list[2])>1:
            flag='indel'
        else:
            flag='snp'
        # b1_index=((b1_geno[p2_num]-(b1_dep/2)-b1_geno[p1_num])/(b1_dep/2))
        # b2_index=((b2_geno[p2_num]-(b2_dep/2)-b2_geno[p1_num])/(b2_dep/2))
        if snp_index==True:
            if b1_dep == 0 or b2_dep == 0: continue
            b1_index=(b1_geno[p1_num])/b1_dep
            b2_index=(b2_geno[p1_num])/b2_dep
            total_index=((b1_geno[p1_num])/b1_dep)-((b2_geno[p1_num])/b2_dep)
            virtual_index=binom_index2(b1_dep,b2_dep,binom_dict)
        elif ED==True:
            if b1_dep == 0 or b2_dep == 0: continue
            b1_index=(b1_geno[p1_num]/b1_dep-b2_geno[p1_num]/b2_dep)**2
            b2_index=(b1_geno[p2_num]/b1_dep-b2_geno[p2_num]/b2_dep)**2
            total_index=b1_index+b2_index
            virtual_index=binom_index(b1_dep,b2_dep,binom_dict)

        # virtual_index=binom_index2(b1_dep/2,b2_dep/2,binom_dict)
        out_snp_list.append(vcf_list+[virtual_index,b1_index,b2_index,total_index,flag])
    return [key,out_snp_list]

def read_dhhp_file(infile):
    loop_file=open(infile,'r')
    out_dict={}
    for line in loop_file:
        line_list=line.strip().split('\t')
        if line_list[0] not in out_dict.keys():
            out_dict[line_list[0]]=[]
        ad_list=[int(line_list[1])]+line_list[2:]
        out_dict[line_list[0]].append(ad_list)
    return out_dict

def cal_dis_func(key,win_size,loop_value):
    step=win_size/10
    win_list=[0,win_size]
    total_list=[]
    start_pos=0
    # this_list=[]
    while win_list[1]<loop_value[-1][0]:
        this_list=[[],[],[],[]]
        for loop_list in loop_value[start_pos:]:
            # print(loop_list)
            i_list=[int(loop_list[0]),float(loop_list[-5]),float(loop_list[-4]),float(loop_list[-3]),float(loop_list[-2])]
            if win_list[0]<=i_list[0]<win_list[1]:
                this_list[0].append(i_list[1])
                this_list[1].append(i_list[2])
                this_list[2].append(i_list[3])
                this_list[3].append(i_list[4])
            elif i_list[0]>=win_list[1]:
                try:
                    total_list.append([(win_list[0]+win_list[1])/2,sum(this_list[0])/len(this_list[0]),sum(this_list[1])/len(this_list[0]),\
                       sum(this_list[2])/len(this_list[0]),sum(this_list[3])/len(this_list[0]), len(this_list[0]),]) 
                except ZeroDivisionError:
                    total_list.append([(win_list[0]+win_list[1])/2,0,0,0,0,0])

                win_list[0]+=step
                win_list[1]+=step
                break
            elif i_list[0]<win_list[0]:
                start_pos+=1
    return [key,total_list]
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='F2_BSA 1.0',\
        epilog='Please cite "OcBSA: an NGS-based Bulk Segregant Analysis Tool for Outcross Populations, Molecular Plant, 2024." Developed by: zhanglk960127@163.com',)
    parser.add_argument('--version', action='version', version='F2_BSA 1.0')
    group = parser.add_mutually_exclusive_group(required=True)
    
    # 添加两个互斥的参数
    group.add_argument('-snpindex','--snpindex', action='store_true', help='Option 1 description')
    group.add_argument('-ED', '--ED', action='store_true', help='Option 2 description')
    parser.add_argument('-p1', '--parent1', type=int, help='Column number of parent1 in the VCF',)
    parser.add_argument('-p2', '--parent2', type=int, help='Column number of arent2 in the VCF',)
    parser.add_argument('-b1', '--pool1', type=int, help='Column number of pool with parent1 trait in the VCF',)
    parser.add_argument('-b2', '--pool2', type=int, help='Column number of pool with parent2 trait in the VCF',)
    parser.add_argument('-d1', '--parentdep1', type=int, help='Minimum coverage of the parents', default=10)
    parser.add_argument('-d2', '--pooldep1', type=int, help='Minimum coverage of the pools度', default=20)
    parser.add_argument('-d3', '--parentdep2', type=int, help='Maximum coverage of the parents度', default=80)
    parser.add_argument('-d4', '--pooldep2', type=int, help='Maximum coverage of the pools', default=200)
    parser.add_argument('-w', '--win', type=int, help='Size of sliding windows',default=1000000)
    parser.add_argument('-vcf', '--input_vcf', help='Path of VCF file ')
    parser.add_argument('-table', '--input_format', help='Path of simple VCF file')
    parser.add_argument('-infile', '--input_infile', help='The intermediate file (ED/snpindex) generated earlier can be used to resize the window; ')
    parser.add_argument('-o', '--out', help='Name of output file; 输出文件名', required=True)
    args = parser.parse_args()
    p1=args.parent1
    p2=args.parent2
    b1=args.pool1
    b2=args.pool2
    snp_index=args.snpindex
    ED=args.ED
    dep_lim=args.parentdep1
    dep_high=args.parentdep2
    dep_lim_pool=args.pooldep1
    dep_high_pool=args.pooldep2
    win_size=args.win
    outfile=args.out    
    infile_vcf=args.input_vcf
    infile_bsa=args.input_infile
    infile_format=args.input_format
    if infile_vcf:
        if p1==None or p2==None or b1==None or b2==None:
            print('You have to provide -p1 -p2 -b1 -b2')
            sys.exit()
    if infile_vcf==None and infile_bsa==None and infile_format==None:
        print('You have to provide .vcf or .bsa file')
        sys.exit()
    ###get binom test dict
    binom_dict=pre_binom_test(10000)
    if infile_bsa:
        print('Start to read the file '+ infile_bsa)
        chr_geno_dict_filter=read_dhhp_file(infile_bsa)
    else:
    # chr_geno_dict={}
        if infile_vcf:
            print('Start reading file '+infile_vcf)
            chr_geno_dict=read_vcf_data(infile_vcf,p1,p2,b1,b2)
        elif infile_format:
            print('Start reading file '+infile_format)
            chr_geno_dict=read_table_file(infile_format)
        chr_geno_dict_filter={}
        ##########filter data
        if snp_index:
            print('Start calculating the snp_index')
        else:
            print('Start calculating the ED value')

        pool = multiprocessing.Pool(len(chr_geno_dict))
        results=[]
        for key,value in chr_geno_dict.items():

            print(key)
            # filter_list=filter_data(key,value,dep_lim,dep_high,dep_lim_pool,dep_high_pool,binom_dict)
            # chr_geno_dict_filter[key]=filter_list
            results.append(pool.apply_async(filter_data, (key,value,dep_lim,dep_high,dep_lim_pool,dep_high_pool,binom_dict,snp_index,ED)))
        pool.close()
        pool.join()
        for result in results:
            chr_list=result.get()
            ########## filter by snp number ###########
            # if len(chr_list[1])<10:
            #     continue
            chr_geno_dict_filter[chr_list[0]]=chr_list[1]
        ####bsa file output
        if snp_index:
            out_bsa_file=open(outfile+'.snpindex','w')
        else:
            out_bsa_file=open(outfile+'.ED','w')

        for key,value in chr_geno_dict_filter.items():
            for i_list in value:
                out_bsa_file.write(key+'\t')
                for i in i_list:
                    out_bsa_file.write(str(i)+'\t')
                out_bsa_file.write('\n')
        out_bsa_file.close()
        if snp_index:
            print('The calculation has been completed and the results are stored in file '+outfile+'.snp')
        else:
            print('The calculation has been completed and the results are stored in file '+outfile+'.ED')

    
    #######window####
    pool = multiprocessing.Pool(len(chr_geno_dict_filter))
    result=[]
    for key,value in chr_geno_dict_filter.items():
        print(key)
        print('Start smoothing  value, window:'+str(win_size))
        result.append(pool.apply_async(cal_dis_func,(key,win_size,value)))

        # lowess_data(key,value)
    pool.close()
    pool.join()
    file3=open(outfile,'w')
    for i_list in result:
        i_list=i_list.get()
        for i in i_list[1][1:]:
            file3.write(i_list[0]+'\t'+str(int(i[0]))+'\t')
            for x in i[1:]:
                file3.write(str(x)+'\t')
            file3.write('\n')
    print('All jobs have been completed and the smoothed results are stored in '+outfile)
