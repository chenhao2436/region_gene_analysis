# -*- coding: utf-8 -*-
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import numpy as np
def Oc_plot(input,out_put,pos,color):
    file1=open(input,'r')
    if pos:
        fig = plt.figure(figsize=(10,5))
    else:
        fig = plt.figure(figsize=(30,5))
    matplotlib.rcParams['pdf.fonttype']=42
    gs = matplotlib.gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0:1,0:1],facecolor='white')

    snp_num_list=[]
    chr_dict={}
    if pos:
        chr=pos.split(',')[0]
        filter_start=int(pos.split(',')[1])
        filter_end=int(pos.split(',')[2])
    for line in file1:
        if line.startswith('#'):
            continue
        line_list=line.strip().split('\t')
        snp_num_list.append(float(line_list[4]))
        # if int(line_list[4])<60 :
        #     continue
        if pos:
            if line_list[0] !=chr:
                continue
            elif int(line_list[1])<filter_start:
                continue
            elif int(line_list[1])>filter_end:
                continue
        
        if line_list[0] in chr_dict.keys():
            chr_dict[line_list[0]].append([float(line_list[2]),float(line_list[3]),float(line_list[4]),int(line_list[1])])
        else:
            chr_dict[line_list[0]]=[[float(line_list[2]),float(line_list[3]),float(line_list[4]),int(line_list[1])]]
    chr_list=sorted(chr_dict.keys(),key=lambda x:int(re.findall('\d+',x)[0]))
    total_list_y=[]
    total_list_y1=[]
    total_list_x=[]
    start=0
    chr_pos_list=[]
    chr_name_list=[]
    chr_pos2_list=[]
    color_list=[]
    big_num=np.percentile(snp_num_list,90)
    for chr in chr_list:
        chr_name_list.append(chr)
        chr_pos2_list.append(start+chr_dict[chr][-1][-1]/2)
        chr_pos_list.append([chr,start+chr_dict[chr][-1][-1]/2,start])
        for i in chr_dict[chr]:
            total_list_x.append(i[3]+start)
            total_list_y.append(i[0])
            total_list_y1.append(i[1])
            if i[2]>big_num:
                color_list.append(big_num)
            else:
                color_list.append(i[2])
        start+=chr_dict[chr][-1][-1]
    for chr in chr_pos_list:
        plt.axvline(chr[2],ls="--",c="black",linewidth='0.3')
    sorted_y_list=sorted(total_list_y1)
    
    # ax1.plot(total_list_x,total_list_y,color='#ed1c24')
    axx=ax1.scatter(total_list_x,total_list_y1,s=8,c=color_list,cmap=color)

    if pos:
        plt.xlim([total_list_x[0],start])
        
        step=round(int((filter_end-filter_start)/20),-6)
        if step==0:
            step=round(int((filter_end-filter_start)/10),-6)
        step=int(step)
        x_ticks = [i for i in range(filter_start,filter_end, step)]
        plt.xticks(x_ticks)
    else:
        top99_values = sorted_y_list[ int(len(sorted_y_list) * 0.99)]
        top95_values = sorted_y_list[ int(len(sorted_y_list) * 0.95)]
        top999_values = sorted_y_list[ int(len(sorted_y_list) * 0.999)]
        ax1.axhline(y=top95_values, color='r', linestyle='--') 
        ax1.axhline(y=top99_values, color='g', linestyle='--') 
        ax1.axhline(y=top999_values, color='b', linestyle='--') 
        plt.xlim([0,start])
        plt.xticks(chr_pos2_list,chr_name_list,fontsize=15)
    plt.yticks(fontsize=15)
    position=fig.add_axes([0.55, 0.01, 0.3, 0.02])
    fig.colorbar(axx,ax=ax1,cax=position,orientation='horizontal')
    plt.savefig(out_put,bbox_inches="tight",pad_inches=0.6,dpi=600)    

def ED_plot(infile,outfile,pos,snp,color):
    file1=open(infile,'r')
    fig = plt.figure(figsize=(30,5))
    matplotlib.rcParams['pdf.fonttype']=42
    gs = matplotlib.gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0:1,0:1],facecolor='white')

    snp_num_list=[]
    chr_dict={}
    if pos:
        chr=pos.split(',')[0]
        filter_start=int(pos.split(',')[1])
        filter_end=int(pos.split(',')[2])
    for line in file1:
        if line.startswith('#'):
            continue
        line_list=line.strip().split('\t')
        snp_num_list.append(float(line_list[-1]))
        # if int(line_list[4])<60 :
        #     continue
        if pos:
            if line_list[0] !=chr:
                continue
            elif int(line_list[1])<filter_start:
                continue
            elif int(line_list[1])>filter_end:
                continue
        # if int(line_list[-1])<60 or 'Un' in line_list[0]:
        #     continue
        if line_list[0] in chr_dict.keys():
            chr_dict[line_list[0]].append([float(line_list[2]),float(line_list[3]),float(line_list[4]),float(line_list[5]),int(line_list[1]),int(line_list[-1])])
        else:
            chr_dict[line_list[0]]=[[float(line_list[2]),float(line_list[3]),float(line_list[4]),float(line_list[5]),int(line_list[1]),int(line_list[-1])]]
    chr_list=sorted(chr_dict.keys(),key=lambda x:int(re.findall('\d+',x)[0]))
    total_list_y=[]
    total_list_y1=[]
    total_list_y2=[]
    total_list_y3=[]
    total_list_x=[]
    start=0
    chr_pos_list=[]
    chr_name_list=[]
    chr_pos2_list=[]
    color_list=[]
    big_num=np.percentile(snp_num_list,90)
    print(big_num)
    for chr in chr_list:
        chr_name_list.append(chr)
        chr_pos2_list.append(start+chr_dict[chr][-1][-2]/2)
        chr_pos_list.append([chr,start+chr_dict[chr][-1][-2]/2,start])
        for i in chr_dict[chr]:
            # print(i)
            total_list_x.append(i[-2]+start)
            total_list_y.append(i[0])
            total_list_y1.append(i[1])
            total_list_y2.append(i[2])
            total_list_y3.append(i[3])
            if i[-1]>big_num:
                color_list.append(big_num)
            else:
                color_list.append(i[-1])
        start+=chr_dict[chr][-1][-2]
    color_list[-1]=0
    for chr in chr_pos_list:
        plt.axvline(chr[2],ls="--",c="black",linewidth='0.3')
    sorted_y_list=sorted(total_list_y1)

    top99_values = sorted_y_list[ int(len(sorted_y_list) * 0.99)]
    top95_values = sorted_y_list[ int(len(sorted_y_list) * 0.95)]
    top999_values = sorted_y_list[ int(len(sorted_y_list) * 0.999)]

    # ax1.plot(total_list_x,total_list_y,color='#ed1c24')
    axx=ax1.scatter(total_list_x,total_list_y3,s=2,c=color_list,cmap=color)
    
    # if snp:
        # ax1.plot(total_list_x,total_list_y1, color='g',lw=0.3, label='Pool1 snp index') 
        # ax1.plot(total_list_x,total_list_y2, color='b', lw=0.3 ,label='Pool2 snp index') 
    # plt.ylim([0.75,0.78])
    
    if pos:
        plt.xlim([total_list_x[0],start])
        step=round(int((filter_end-filter_start)/20),-6)
        if step==0:
            step=round(int((filter_end-filter_start)/10),-6)
        x_ticks = [i for i in range(filter_start,filter_end, step)]
        plt.xticks(x_ticks)
    else:
        ax1.plot(total_list_x,total_list_y, color='r',lw=0.5 ,label='thresholds') 
        plt.xlim([0,start])
        plt.xticks(chr_pos2_list,chr_name_list,fontsize=15)
    ax1.legend()
    ax1.legend(loc='upper right')
    plt.yticks(fontsize=15)
    position=fig.add_axes([0.55, 0.01, 0.3, 0.02])
    fig.colorbar(axx,ax=ax1,cax=position,orientation='horizontal')
    plt.savefig(outfile,bbox_inches="tight",pad_inches=0.6,dpi=600)    



if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='BSA-fig 1.0',\
        epilog='Please cite "OcBSA: an NGS-based Bulk Segregant Analysis Tool for Outcross Populations, Molecular Plant, 2024." Developed by: zhanglk960127@163.com',)
    parser.add_argument('--version', action='version', version='BSA-fig 1.0')
    parser.add_argument('-f', '--inputfile', type=str, help='input file for plot, output_file of OcBSA or F2_BSA; ',required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-OcValue','--OcValue', action='store_true', help='plot for OcValue', )
    
    group.add_argument('-snpindex','--snpindex', action='store_true', help='plot for snpindex')
    group.add_argument('-ED', '--ED', action='store_true', help='plot for ED')
    parser.add_argument('-p','--position', type=str, help='Select a coordinate to plot a portion of it. e.g. chr10,1,10000')
    parser.add_argument('-c','--color', type=str, help='Choose a gradient color (color of heatmap) for the dot plot', default='plasma_r')
    parser.add_argument('-o', '--output', type=str, help='The name of the output figure, ending with .png or .pdf ',required=True)
    args = parser.parse_args()
    infile=args.inputfile
    OcValue=args.OcValue
    snpindex=args.snpindex
    ED=args.ED
    outfile=args.output
    color=args.color
    pos=args.position
    if outfile.endswith('.png') or outfile.endswith('.pdf') :
        if OcValue:
            Oc_plot(infile,outfile,pos,color)
        elif ED:
            ED_plot(infile,outfile,pos,snpindex,color)
        elif snpindex:
            ED_plot(infile,outfile,pos,snpindex,color)






    


   
