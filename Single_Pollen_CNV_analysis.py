import os,glob,math,sys
import pysam
from multiprocessing import Pool


sample_list = ['CAU5P10S1', 'CAU5P10S2', 'CAU5P1S1', 'CAU5P1S2', 'CAU5P3S1', 'CAU5P3S2', 'CAU5P6S1', 'CAU5P6S2', 'CAU5P8S1', 'CAU5P8S2']
bam_list = ['CAU5P10S1.sort.bam', 'CAU5P10S2.sort.bam', 'CAU5P1S1.sort.bam', 'CAU5P1S2.sort.bam', 'CAU5P3S1.sort.bam', 'CAU5P3S2.sort.bam', 'CAU5P6S1.sort.bam', 'CAU5P6S2.sort.bam', 'CAU5P8S1.sort.bam', 'CAU5P8S2.sort.bam']

#sample_list = ['CAU5P1S1']
#bam_list = ['CAU5P1S1.sort.bam']

sbin_list = []
step_bin_size = 5000
bin_size = 500000

chrs = [['1',307041717],['2',244442276],['3',235667834],['4',246994605],['5',223902240],['6',174033170],['7',182381542],['8',181122637],['9',159769782],['10',150982314]]

def main():
    print 'step1'
    sbin_list = [[] for each in sample_list]
    p_out = []
    p = Pool(processes=10)
    for i in range(len(sample_list)):
        p_out.append([])
        p_out[i] = p.apply_async(sample_depth, args=(sample_list[i], bam_list[i], step_bin_size, chrs,))
    p.close()
    p.join()
    sbin_list = []
    for each in p_out:
        sbin_list.append(each.get())
#    print sbin_list
    for sp in sbin_list:
        print '======='
        for chr in sp:
            print len(chr)

    print 'step2'
    refined_bin = refine_bin_size(sbin_list, step_bin_size, bin_size)
    print 'step3'
    total_CNR_list = []
    for i in range(len(sample_list)):
        CNR_list = cal_CNR(sample_list[i], sbin_list[i], step_bin_size, refined_bin)
        o = open(sample_list[i]+'.xls','w')
        for each in CNR_list:
             o.write('\t'.join(map(str,each))+'\n')
        o.close()
        total_CNR_list.append(CNR_list)
    print 'step4'
    plot_CNR(total_CNR_list,sample_list)

def sample_depth(sp, bamfile, step_bin_size, chrs):
    print sp
    inbam = pysam.AlignmentFile(bamfile, "rb")
    sbin = []
    for chr,chr_l in chrs:
        print sp,chr
        sbin.append([0]*(chr_l/step_bin_size+1))
        i = 0
        for read in inbam.fetch(chr, 0, chr_l):
            while step_bin_size*(i+1) < read.reference_start:
                i += 1
            sbin[-1][i] += 1
    return sbin

def refine_bin_size(sbin_list, step_bin_size, bin_size):
    total_list = []
    total_bin_reads = 0
    total_bin_number = 0
    for i in range(10):
        total_list.append([])
        big_index = -1
        for j in range(len(sbin_list[0][i])):
            if j/100 > big_index:
                big_index += 1
                total_list[-1].append(0)
            for eachsp in sbin_list:
                total_list[-1][-1] += eachsp[i][j]
        total_bin_reads += sum(total_list[-1][:-1])
        total_bin_number += (len(total_list[-1])-1)
#    print total_list

    average_bin_depth = total_bin_reads/total_bin_number
    print 'average_bin_depth',average_bin_depth
    refined_bin = []
    for i in range(10):
        refined_bin.append([0])
        tmp_depth = 0
        for j in range(len(sbin_list[0][i])):
            for eachsp in sbin_list:
                tmp_depth += eachsp[i][j]
            refined_bin[-1][-1] += step_bin_size
            if tmp_depth >= average_bin_depth:
                refined_bin[-1].append(0)
                tmp_depth = 0
#    print 'refined_bin'
#    print refined_bin
    return refined_bin

def cal_CNR(sample, sbin, step_bin_size, refined_bin):
    refined_count = []
    all_for_median = []
    for i in range(10):
        refined_count.append([])
        g = (x for x in sbin[i])
        for each_bin in refined_bin[i]:
            l = each_bin/step_bin_size
            refined_count[-1].append(0)
            for j in range(l):
                refined_count[-1][-1] += g.next()
        all_for_median += refined_count[-1][:-1]
    median_count = get_median(all_for_median)
    CNR_list = []
    for i in range(10):
        CNR_list.append([])
        CNR_list[i] = [math.log((float(x)/median_count),2) for x in refined_count[i][:-1]]
#        CNR_list[i] = [math.log(x/median_count, 2) for x in refined_count[i]]
    return CNR_list

def get_median(data):
   data.sort()
   half = len(data)//2
   return (data[half] + data[~half])/2

def plot_CNR(total_CNR_list,sample_list):
    colors = ['rgb(30,144,255)','rgb(255,0,255)']
    txt = '''
    <svg version="1.1"
    baseProfile="full"
    width="2000" height="3000"
    xmlns="http://www.w3.org/2000/svg">
'''
    for i in range(len(total_CNR_list)):
        CNR_list = total_CNR_list[i]
        x_stand = 200
        y_stand = 50+i*200
        txt += '''
    <text x="{xsp}" y="{ysp}" font-size="16" text-anchor="middle" transform="rotate(-90 {xsp},{ysp})" fill="black">{sp}</text>
    <text x="{xt}" y="{ysp}" font-size="14" text-anchor="middle" transform="rotate(-90 {xt},{ysp})" fill="black">log2(NCR)</text>
    <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" style="stroke:rgb(0,0,0);stroke-width:2;stroke-opacity:0.5"/>
    <line x1="{x3}" y1="{y3}" x2="{x3}" y2="{y4}" style="stroke:rgb(0,0,0);stroke-width:2;stroke-opacity:0.5"/>
    <line x1="{x4}" y1="{y8}" x2="{x5}" y2="{y8}" style="stroke:rgb(0,0,0);stroke-width:2;stroke-opacity:0.5"/>
    <text x="{x8}" y="{y5}" text-anchor="middle" fill="black">2</text>
    <line x1="{x6}" y1="{y9}" x2="{x7}" y2="{y9}" style="stroke:rgb(0,0,0);stroke-width:2;stroke-opacity:0.5"/>
    <text x="{x9}" y="{y6}" text-anchor="middle" fill="black">-2</text>
    <line x1="{x6}" y1="{y1}" x2="{x7}" y2="{y1}" style="stroke:rgb(0,0,0);stroke-width:2;stroke-opacity:0.5"/>
'''.format(xsp=80,ysp=y_stand+100, sp=sample_list[i], xt=120, x1=x_stand, y1=y_stand+100, x2=x_stand+820, y2=y_stand+100, x3=180, y3=y_stand+100-36, y4=y_stand+100+36, x4=170, x5=180, y5=y_stand+100-36+4, x8=160, x6=170, x7=180, y6=y_stand+100+36+4, x9=160, y8=y_stand+100-36, y9=y_stand+100+36)
        c = 0
        j = 0
        for chr_list in CNR_list:
            color = colors[c%2]
            c += 1
            for p in chr_list:
                if abs(p) < 4:
                    x = j/5.0 + x_stand
                    y = -p*18 + y_stand + 100
                    txt += '''
    <circle cx="{x}" cy="{y}" r="1" fill="{color}"/>
'''.format(x=x, y=y, color=color)
                j += 1   
    open('pollen_deletion.svg','w').write(txt+'\n    </svg>\n')

if __name__ == '__main__':
    main()


