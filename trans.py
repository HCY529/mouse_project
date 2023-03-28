def run(dirs_path):
    f = open(dirs_path , 'r')
    name = dirs_path.split('/')[-1].split('.')[0]
    w = open(name+'.dss.input.txt','w')
    w.write('chr'+'\t'+'pos'+'\t'+'N'+'\t'+'X'+'\n')
    for line in f:
        d = line.strip().split('\t')
        col_3 = int(d[-2])+int(d[-1])
        w.write(d[0]+'\t'+d[1]+'\t'+str(col_3)+'\t'+d[-1]+'\n')
    f.close()
    w.close()

run("/home/huchenyuan/mouse_project/data/methylation_information_2/SRR833527_1_val_1_bismark_bt2_pe.bismark.cov")
run("/home/huchenyuan/mouse_project/data/methylation_information_2/SRR833528_1_val_1_bismark_bt2_pe.bismark.cov")
