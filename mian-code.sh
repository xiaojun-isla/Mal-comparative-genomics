############################################松墨天牛基因组比较基因组学分析分析全代码##################################################
###整合后的代码按照逻辑顺序来书写，具体代码内的路径不一致按照分析逻辑进行修改。

1.提取不含isoform的物种蛋白质序列文件。软件：TBtools、
##先提取第三列为CDS的行的第九列
awk -F "\t" '$3 == "CDS" {print $9}' genomic.gff > output.txt

##对于上一步得到的output文件删除掉一模一样的行
sort output.txt | uniq > output_unique.txt
##成功得到output_unique.txt文件

##对于id.txt文件中的ID，每次到output_unique.txt中查找，输出匹配行，对于匹配行，输出protein_id开始，到；为止的字符串
dos2unix id.txt
for id in $(cat id.txt);
do
##echo $id
grep $id output_unique.txt | grep -o 'protein_id=[^;]*' >>final_id.txt  #grep $id output_unique.txt | grep -o 'ID=cds.[^;]*' >>final_id.txt
done
sort final_id.txt | uniq > final_ids.txt





2.基因家族聚类。软件：Orthofinder
#!/bin/bash
#PBS -N ortho
#PBS -l nodes=6:ppn=3
#PBS -q blade
#PBS -l walltime=999:00:00
#PBS -j oe

##source /miniconda3/envs/ortho/bin/orthofinder
/miniconda3/envs/ortho/bin/orthofinder -f /malgenome/zfiles/Lepidoptera/2ortho/faafile/ -M msa -t 16 -a 16 -S diamond





3.基因家族的收缩和扩张分析
#!/bin/bash
#PBS -N cafe
#PBS -l nodes=6:ppn=3
#PBS -q blade
#PBS -l walltime=999:00:00
#PBS -j oe

##先把cafe相关的脚本下载下来
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafecore.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_clade_and_size_filter.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_draw_tree.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_longest_iso.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_mcl2rawcafe.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_prep_r8s.py
wget https://github.com/hahnlab/cafe_tutorial/blob/main/python_scripts/cafetutorial_report_analysis.py

##准备r8s的输入文件r8s_ctl_file.txt
##seqkit可以提取物种树比对的氨基酸序列数目（？这里得到的数字应该是用于推断物种数的比对序列的碱基数目）
seqkit stat /malgenome/zfiles/Lepidoptera/3cafe/SpeciesTreeAlignment.fa
python2 /software/CAFE5/tutorial/prep_r8s.py -i /malgenome/zfiles/Lepidoptera/3cafe/SpeciesTree_rooted.txt -o /malgenome/zfiles/Lepidoptera/3cafe/r8s_ctl_file.txt -s 985613 -p 'Diabrotica_virgifera,Photinus_pyralis' -c '262'

##运行r8s，计算超度量树
/software/r8s/r8s1.81/src/r8s -b -f /malgenome/zfiles/Lepidoptera/3cafe/r8s_ctl_file.txt >r8s_tmp.txt
tail -n 1 /malgenome/zfiles/Lepidoptera/3cafe/r8s_tmp.txt | cut -c 16- > /malgenome/zfiles/Lepidoptera/3cafe/r8s_ultrametric.txt

##把Orthogroups.GeneCount.tsv改成cafe需要的文件格式
awk 'OFS="\t" {$NF=""; print}' /malgenome/zfiles/5cafe/Orthogroups.GeneCount.tsv > tmp && awk '{print "(null)""\t"$0}' tmp > /malgenome/zfiles/5cafe/cafe.input.tsv && sed -i '1s/(null)/Desc/g' /malgenome/zfiles/5cafe/cafe.input.tsv && rm tmp

##去除基因拷贝变异数特别大的基因家族。该脚本可以过滤掉一个或多个物种有超过100个基因拷贝的基因家族
python2 /software/CAFE5/tutorial/clade_and_size_filter.py -i cafe.input.tsv -o filtered.cafe.input.tsv -s
##这里运行的时候不需要加全路径，不然会出错

##文件都已经准备就绪，可以开始跑了
#/miniconda3/envs/mal/bin/cafe5 -i /malgenome/zfiles/5cafe/filtered.cafe.input.tsv -t /malgenome/zfiles/5cafe/r8s_ultrametric.txt -p -k 2 -o /malgenome/zfiles/5cafe/cafe5 -c 10 
/software/CAFE5/bin/cafe5 -i /malgenome/zfiles/5cafe/filtered.cafe.input.tsv -t /malgenome/zfiles/5cafe/r8s_ultrametric.txt -p -k 2 -o /malgenome/zfiles/5cafe/cafe5.1 -c 18
#如果报了类似Tribolium_castaneum was not found in gene family OG0000000的错，用dos2unix修复一下Hiefiltered.cafe.input.tsv文件

#cafe结果文件说明。/malgenome/zfiles/Lepidoptera/3cafe/
Gamma_asr.tre 每个基因家族的树文件，树上的星号代表这个物种的基因家族有显著变化
Gamma_branch_probabilities.tab 每个分支计算的概率
Gamma_category_likelihoods.txt 
Gamma_change.tab 每个基因家族在每个节点的收缩扩张数目
Gamma_clade_results.txt 进化树上每个节点扩张或者收缩了多少基因家族
Gamma_count.tab 每个基因家族在每个节点的基因数
Gamma_family_likelihoods.txt 
Gamma_family_results.txt 基因家族变化的p值和是否显著的结果，0.05
Gamma_report.cafe cafe5没有这个文件，导致不能直接作图，CAFE_fig.py没有了输入文件，因此目前还未正式release的cafe5.1经过设置后又能输出这个结果文件了。
Gamma_results.txt 模型，最终似然值和最终lambda值等参数信息





4.系统发育分析
#!/bin/bash
#PBS -N raxml
#PBS -l nodes=1:ppn=10
#PBS -q blade
#PBS -l walltime=9999:00:00
#PBS -j oe

##物种树。Orthofinder中输出的那个SpeciesTree_rooted.txt并不是物种树，不是用单拷贝同源基因序列构建的，因此要用单拷贝同源基因再做一个物种发育分析。
#在gxx环境下工作
/miniconda3/envs/gxx/bin/trimal -in SpeciesTreeAlignment.fa -out SpeciesTreeAlignment_trimal.fa -fasta -gt 0.6 -cons 60
cd /malgenome/zfiles/Lepidoptera/4phylogeny/
/miniconda3/envs/gxx/bin/raxmlHPC-PTHREADS -T 20 -m PROTGAMMAJTT -f a -p 123 -x 123 -# 100 -k -o Drosophila_melanogaster -n out -s /malgenome/zfiles/Lepidoptera/4phylogeny/SpeciesTreeAlignment_trimal.fa
##准备r8s的输入文件r8s_ctl_file.txt
##seqkit可以提取物种树比对的氨基酸序列数目
seqkit stat /malgenome/zfiles/Lepidoptera/4phylogeny/SpeciesTreeAlignment.fa
python2 /software/CAFE5/tutorial/prep_r8s.py -i /malgenome/zfiles/Lepidoptera/4phylogeny/RAxML_bestTree.out -o /malgenome/zfiles/Lepidoptera/4phylogeny/r8s_ctl_file.txt -s 985613 -p 'Diabrotica_virgifera,Photinus_pyralis' -c '262'
##运行r8s，计算超度量树
/software/r8s/r8s1.81/src/r8s -b -f /malgenome/zfiles/Lepidoptera/4phylogeny/r8s_ctl_file.txt >r8s_tmp.txt
tail -n 1 /malgenome/zfiles/Lepidoptera/4phylogeny/r8s_tmp.txt | cut -c 16- > /malgenome/zfiles/Lepidoptera/4phylogeny/r8s_ultrametric.txt





5.KEGG & GO分析
#(1)节点基因的提取
#提取显著扩张或收缩的 orthogroups ID
cat Gamma_family_results.txt |grep "y"|cut -f1 >p0.05.significant

# 显著扩张 / 收缩的基因家族在每个节点的收缩与扩增数目
grep -f p0.05.significant Gamma_change.tab > Gamma_p0.05change.tab

#提取Gamma_p0.05change.tab第14列松墨mal的显著扩张的orthogroupsID
awk '$2 ~ /^[1-9][0-9]*$/ {print $1}' Gamma_p0.05change.tab > agr.expanded  正整数前面没有+

#提取Gamma_change.tab第8列松墨mal的收缩的orthogroupsID
#cat Gamma_change.tab |cut -f1,8 |grep "-" >mal.contracted

#提取显著扩张的基因列表。
#松墨天牛
grep -f mal.expanded Orthogroups.txt |sed "s/ /\n/g"|grep "evm.model.Contig" |sort -k 1.17n |uniq >mal.0.05.expanded.genes 
#其他叶节点
awk 'NR==FNR{ids[$1]; next} $1 in ids' agl.expanded Orthogroups.tsv > agl_tmp.txt #提取含有og的行，用Excel打开agl_tmp.txt，删除多余列和第一列，然后把, 替换成，
awk '{gsub(/,/, "\n"); print > "agl.genes1"}' agl_tmp.txt #将每行内容按,拆分成多行
sed 's/"//g' agl.genes1 | awk '{gsub(/,/, "\n"); print > "agl.genes2"}' #删除多余的"，然后用Excel打开agl.genes2，删掉除了ID之外的多余内容，Excel里用分列功能
sort agl.genes2 | uniq > agl.0.05.expanded.genes
最后删除：agl_tmp.txt agl.genes1 agl.genes2
#其他节点
awk 'NR==FNR{ids[$1]; next} $1 in ids' node30.expanded Orthogroups.tsv > node30_tmp.txt #提取含有og的行，用Excel打开agl_tmp.txt，删除多余列和第一列，只保留该节点所有子节点的列以及第一列，然后把, 替换成，。然后将所有内容放到一列里
awk '{gsub(/,/, "\n"); print > "node30.genes1"}' node30_tmp.txt #将每行内容按,拆分成多行
sed 's/"//g' node30.genes1 | awk '{gsub(/,/, "\n"); print > "node30.genes2"}' #删除多余的"，然后用Excel打开agl.genes2，删掉除了ID之外的多余内容，Excel里用分列功能
sort node30.genes2 | uniq > node30.0.05.expanded.genes
最后删除：node24_tmp.txt node24.genes1 node24.genes2

#提取显著收缩的基因列表。########################收缩暂时有争议未解决，但是可以先提取
#grep -f mal.0.05.contracted Orthogroups.txt | sed "s/ /\n/g"|grep "evm.model.Contig" |sort -k 1.17n |uniq >mal.0.05.contracted.genes 如果无法成功可以结合Excel手动操作
#grep -f mal.0.01.contracted Orthogroups.txt | sed "s/ /\n/g"|grep "evm.model.Contig" |sort -k 1.17n |uniq >mal.0.01.contracted.genes

#提取显著扩张的蛋白序列
#seqkit grep -f mal.0.05.expanded.genes Monochamus_alternatus.faa >mal.0.05.expanded.pep.fa 

#提取显著收缩的蛋白序列
#seqkit grep -f mal.0.05.contracted.genes Monochamus_alternatus.faa >mal.0.05.contracted.pep.fa

#(2)背景基因的选取。选取松墨天牛节点上的所有基因
#先对pep文件序列名称行进行修改
#sed 's/ .*$//' Anoplophora_glabripennis.faa > modified.faa

#松墨节点中存在基因的 orthogroupsID
awk 'NR!=1 && $9>0 {print $0}' Gamma_count.tab | cut -f1 > mal.ogs


#提取这些OG的基因
#松墨天牛
grep -f mal.ogs Orthogroups.txt |sed "s/ /\n/g" | grep "evm.model.Contig" | sort | uniq > mal.ogs.genes
#其他叶节点
awk 'NR==FNR{ids[$1]; next} $1 in ids' agl.ogs Orthogroups.tsv > agl_output.txt #提取含有og的行，用Excel打开agl_tmp.txt，删除多余列和第一列，然后把, 替换成，
awk '{gsub(/,/, "\n"); print > "agl.ogs.genes1"}' agl_output.txt #将每行内容按,拆分成多行
sed 's/"//g' agl.ogs.genes1 | awk '{gsub(/,/, "\n"); print > "agl.ogs.genes2"}' #删除多余的"，然后用Excel打开rma.ogs.genes2，删掉除了ID之外的多余内容，Excel里用分列功能
sort agl.ogs.genes2 | uniq > agl.ogs.genes
最后删除： agl_output.txt agl.ogs.genes1 agl.ogs.genes2
#其他节点
awk 'NR==FNR{ids[$1]; next} $1 in ids' node30.ogs Orthogroups.tsv > node30_output.txt #提取含有og的行，用Excel打开agl_output.txt，删除多余列和第一列，只保留该节点所有子节点的列以及第一列，然后把, 替换成，。然后将所有内容放到一列里
awk '{gsub(/,/, "\n"); print > "node30.ogs.genes1"}' node30_output.txt #将每行内容按,拆分成多行
sed 's/"//g' node30.ogs.genes1 | awk '{gsub(/,/, "\n"); print > "node30.ogs.genes2"}' #删除多余的"，然后用Excel打开rma.ogs.genes2，删掉除了ID之外的多余内容，Excel里用分列功能
sort node30.ogs.genes2 | uniq > node30.ogs.genes
最后删除： node24_output.txt node24.ogs.genes1 node24.ogs.genes2


#提取序列
seqkit grep -f mal.ogs.genes Monochamus_alternatus.faa >mal.ogs.pep.fa

#用网页版eggnog-mapper对背景基因集注释
mal.ogs.pep.fa

#获取GO的注释信息并修饰
grep "^id:" go-basic.obo |awk '{print $2}' > GO.id
grep "^name:" go-basic.obo |awk '{$1=""; print}' > GO.name 
grep "^namespace:" go-basic.obo |awk '{print $2}' > GO.class
paste GO.id GO.name GO.class -d "\t" > GO.library





6.重复序列和转座子分析
#!/bin/bash
#PBS -N masker
#PBS -l nodes=4:ppn=4
#PBS -q blade
#PBS -l walltime=999:00:00
#PBS -j oe

#(1).repeatmasker+repeatmodeler
#由于直接用repeatmasker的结果不好，重复序列占比太低，甚至不足10%，因此加上repeatmodeler
#1.把碱基全部改成大写的
seqkit seq -u lower.geno >upper.geno
#seqkit seq -u Aethina_tumida.fna >Atu.fna
#seqkit seq -u Dme.fna > Dme.fasta
##确定genome文件中是否有小写字母
#grep '[a-z]' genome.faa

#(2).RepeatModeler
#(a)导出同源物种重复序列库
#在/software/RepeatMasker/下，tree环境下运行
#查找目标物种的最近分支，由于所有物种都是鞘翅目下的，所以使用Coleoptera作为最近分支，这样可以重复使用
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 lineage -ad Coleoptera
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 lineage -ad Lepidoptera
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 lineage -ad Diptera
# 查找近缘物种及其上祖先节点，其下所有类群repeat famlies，输出格式embl。 -a ancestor，-d descendent
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 families -f embl -a -d Coleoptera > Coleoptera.embl	
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 families -f embl -a -d Lepidoptera > Lepidoptera.embl
python famdb.py -i /software/RepeatMasker/Libraries/RepeatMaskerLib.h5 families -f embl -a -d Diptera > Diptera.embl
# 转换格式为fasta，方便后续合并
buildRMLibFromEMBL.pl Coleoptera.embl > Coleoptera.fasta	
buildRMLibFromEMBL.pl Lepidoptera.embl > Lepidoptera.fasta
buildRMLibFromEMBL.pl Diptera.embl > Diptera.fasta

#此时进入/software/RepeatModeler-2.0.4/目录下，可以在te环境下运行，也可以不，记得先建好一个目录，用来存放RepeatModeler的预测结果
# 用基因组组装结果构建数据库
BuildDatabase -name /software/RepeatModeler-2.0.4/Ead/Ead -engine ncbi /malgenome/zfiles/Lepidoptera/6te/db/Ead.fasta 

#(c)self-training （重装所有软件之后脚本提交运行），先安装repeatmasker1.0.5，再安装repeatmodeler2.0.4。这里是在用repeatmodeler de novo预测物种的重复序列家族。
RepeatModeler -threads 16 -database /software/RepeatModeler-2.0.4/Ead/Ead -engine ncbi -quick 1003 结束
#(d)整合数据库。(a)得到的同源物种数据库库 + (c)得到的训练结果sample.families.fa
cat Coleoptera.fasta Dva-families.fa > final_Dva_repeat.fasta
#(3)RepeatMasker
#source /miniconda3/etc/profile.d/conda.sh
#conda activate masker
#cd /software/RepeatMasker/
#/malgenome/zfiles/Lepidoptera/6te/db/
#重装所有软件之后，可以提交脚本运行了
RepeatMasker -nolow -no_is -norna -a -inv -html -gff -engine ncbi -parallel 8 -lib final_Ead_repeat.fasta /malgenome/zfiles/Lepidoptera/6te/db/Ead.fasta -dir /malgenome/zfiles/Lepidoptera/6te/Ead 结束

#(4)计算repeat landscape
./calcDivergenceFromAlign.pl -s Agl.divsum Agl.align
./createRepeatLandscape.pl -div example.divsum > /home/user/public_html/example.html
#在/malgenome/zfiles/Lepidoptera/6te/下提交脚本
/malgenome/zfiles/Lepidoptera/6te/calcDivergenceFromAlign.pl -s Agl.divsum Agl.fasta.align

#在各自的目录下直接提交命令
#可以这样画图
/malgenome/zfiles/Lepidoptera/6te/createRepeatLandscape.pl -g 706968555 -div Agl.divsum > Agl.html
#也可以直接从divsum文件提取，然后用Excel修改表格后用R画图
tail -n 72 Agl.divsum > Agl.Kimura.distance





7.结构域分析
#!/bin/bash
#PBS -N pfam-scan
#PBS -l nodes=1:ppn=10
#PBS -q fat
#PBS -l walltime=999:00:00
#PBS -j oe

#下载pfam-A的网址
#http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/

#(1).下载库
wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget ftp://ftp.ebi.ac.uk:21/pub/databases/Pfam/current_release/active_site.dat.gz
gunzip *.gz

#(2).建库
hmmpress Pfam-A.hmm

#(3).使用pfam-scan进行基因家族鉴定，初步确定OG的功能。但是要先对cafe的结果进行统计，先统计出感兴趣的节点或者物种收缩和扩张的基因家族，然后才能做KG或者pfam
#此时在/malgenome/zfiles/2pfam_scan/目录下，在pfam环境下工作
source /miniconda3/bin/activate /miniconda3/envs/pfam

for file in mal.0.01.contracted.pep.fa mal.0.01.expanded.pep.fa mal.0.05.contracted.pep.fa mal.0.05.expanded.pep.fa
do
/miniconda3/envs/pfam/bin/pfam_scan.pl \
-fasta /malgenome/zfiles/2pfam_scan/mal-seq/${file} \
-dir /malgenome/zfiles/2pfam_scan/pfam-A \
-outfile /malgenome/zfiles/2pfam_scan/${file}.txt
done





8.基因组共线性分析。软件：使用过lastz、MCScanX、MCScan。要在jcvi环境中运行
#!/bin/bash
#PBS -N mcscanx
#PBS -l nodes=5:ppn=3
#PBS -q blade
#PBS -l walltime=9999:00:00
#PBS -j oe

#MCScan
#在做不同物种间的比较时，要在各自的文件夹下工作，并且要在jcvi环境下
#(1).将gff压缩文件直接转换成bed格式（解不解压结果都一样。根据物种的不同，有的时候type可以选择gene）
python -m jcvi.formats.gff bed --type=mRNA --key=ID Monochamus_alternatus.gff3 -o Ma.bed
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id --primary_only Rutpela_maculata.gff3 -o Rm.bed
python -m jcvi.formats.gff bed --type=mRNA --key=Dbxref --primary_only Anoplophora_glabripennis.gff -o Ag.bed
python -m jcvi.formats.gff bed --type=mRNA --key=locus_tag Rhamnusium_bicolor.gff -o Rb.bed
python -m jcvi.formats.gff bed --type=mRNA --key=locus_tag Exocentrus_adspersus.gff -o Ea.bed

#(2).对cds文件进行reformat
python -m jcvi.formats.fasta format Rutpela_maculata.cds.fa Rm.cds
python -m jcvi.formats.fasta format Tianniu_longest.cds.seq.fasta Ma.cds
python -m jcvi.formats.fasta format Ag.cds.fna Ag.cds
python -m jcvi.formats.fasta format cds_from_genomic.fna Ea.cds
python -m jcvi.formats.fasta format cds_from_genomic.fna Rb.cds
sed '/^>/ s/\(.*\)\.1.*/\1/' Rm.cds > new_Rm.cds

#(3).整理bed文件和cds文件，这两个文件中的序列ID要一致、数目也要一致
sed -E '/^>/s/>[^>]*locus_tag=/>/g' Rb.cds > Rb1.cds #删除>和locus_tag=之间的内容
sed -E '/^>/s/\].*$//' Rb.cds > Rb1.cds #删除]及其之后的内容

#(4).开始分析
python -m jcvi.compara.catalog ortholog Ma Rm --no_strip_names --cscore=.99  #共线性分析，结果文件是Ma.Rm.pdf，一个点图，查看两个基因组之间的共线性关系
python -m jcvi.compara.catalog ortholog Ag Rm --no_strip_names --cscore=.99
python -m jcvi.compara.catalog ortholog Ea Rm --no_strip_names --cscore=.99
python -m jcvi.compara.catalog ortholog Rb Rm --no_strip_names --cscore=.99
python -m jcvi.graphics.dotplot Ma.Rm.anchors  #点图，查看成对同线性，和上一条命令生成的图片一样
python -m jcvi.compara.synteny depth --histogram Ma.Rm.anchors  #可以检测用来测试的同线性模型是否是1:1，结果文件Ma.Rm.depth.pdf
python -m jcvi.compara.synteny depth --histogram Ag.Rm.anchors
python -m jcvi.compara.synteny depth --histogram Ea.Rm.anchors
python -m jcvi.compara.synteny depth --histogram Rb.Rm.anchors

#(5).染色体水平的局部同线性图绘制（但是这种图只适合比对双方都是染色体水平的分析，现在的contig数据做不了）
#(a)先生成simple文件
python -m jcvi.compara.synteny screen --minspan=30 --simple Ma.Rm.anchors Ma.Rm.anchors.new  
#(b)准备好seqid和layout文件
#(c)运行绘图模块
python -m jcvi.graphics.karyotype seqid layout  #共线性图

#(6).微共线性可视化（可以用scaffold或contig level）
#(a)得到block文件
python -m jcvi.compara.synteny mcscan Rm.bed Ma.Rm.lifted.anchors --iter=1 -o Ma.Rm.i1.blocks
#(b)筛选要画出来的block文件，这里只画了斑点X染色体的比对结果，将Xblocks整理成要用的形式
sed -n '32101,33598p' Rb.Rm.i1.blocks > Xblocks
awk '{print $2}' Xblocks > Xblocks-2
sort Xblocks-2 | uniq > Xblocks-22
#(c)整理bed文件
python -m jcvi.formats.bed merge Rm.bed Ma.bed -o Rm_Ma.bed
#(d)画图 --format svg
python -m jcvi.graphics.synteny Xblocks Rm_Ma.bed blocks.layout 
##blocks.layout例子如下
# x,   y, rotation,   ha,     va,   color, ratio,            label
0.5, 0.6,        0, left, center,       m,     1,       Rm Chr1
0.5, 0.4,        0, left, center, #fc8d62,     1, peach scaffold_1
# edges
e, 0, 1





9.基因家族分析。
*****提取所有物种基因中的P450基因。*****
#(1)先用blastp提取符合evalue cutoff的物种内的P450基因。
#!/bin/bash
#PBS -N mal
#PBS -l nodes=6:ppn=3
#PBS -q blade
#PBS -l walltime=9999:00:00
#PBS -j oe

#(a)对物种蛋白质文件去重，在mal环境下运行
seqkit rmdup mal.faa -s -i -o mal.clean.faa

#(b)去重之后开始blastp
makeblastdb -in mal.clean.faa -dbtype prot -out mal -parse_seqids
blastp -db /malgenome/zfiles/Lepidoptera/9genefamily/P450/1blastp/maldb/mal -query /malgenome/zfiles/Lepidoptera/9genefamily/P450/all.p450.pep.fa -num_threads 18 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/P450/mal.blastp.out 结束

#(2)pfam结构域
#以基因家族的种子文件作为索引，获得hmm文件
hmmbuild A.hmm PF00067_seed.txt 

#在某个物种中，以上一步获取得到的hmm文件为索引，搜寻基因家族成员，以A.out的结果输出
hmmsearch --domtblout mal.hmm.out A.hmm mal.clean.faa

#(3)blastp和hmm的结果取交集或者取合集，得到每个物种的结果之后，提取相应的氨基酸序列
seqkit grep -f agl agl.clean.faa > agl.p450.pep

#(4)先做MSA，此时在tree环境下工作
#!/bin/bash
#PBS -N muscle
#PBS -l nodes=5:ppn=4
#PBS -q blade
#PBS -l walltime=9999:00:00
#PBS -j oe

~/miniconda3/envs/tree/bin/muscle -align agl.mal.tca.fasta -output agl.mal.tca.align.fasta

#(5)用多序列比对结果建树
#先修改一下序列的ID
awk 'FNR==NR{a[">"$1]=$2;next} /^>/{print ">"a[$0];next} 1' mal mal.p450.pep > mal.p450.pep1 #将mal的p450蛋白序列文件中的名称换成自命名的名称

sed '/^>/ s/ .*$//' agl.p450.pep > agl.p450.pep1 #将agl和tca蛋白序列文件中名称多余的东西删除，只留下accession number，方便修改成自命名的名称
awk 'FNR==NR{a[">"$1]=$2;next} /^>/{print ">"a[$0];next} 1' agl agl.p450.pep1 > agl.p450.pep

#开始建树
~/miniconda3/envs/tree/bin/iqtree2 -s agl.mal.tca.align.afa -nt 20


*****提取所有物种中的PF00089*****
#(1)blastp
#去重的clean文件和建好的库db可以直接从P450下的目录里复制
#直接开始blastp
blastp -db /malgenome/zfiles/Lepidoptera/9genefamily/trypsin/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/trypsin/mal.pf00089.pep.fa -num_threads 18 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/trypsin/1blastp/out/agl.blastp.out 

#(2)pfam结构域
#以基因家族的种子文件作为索引，获得hmm文件
hmmbuild B.hmm SEED.ann

#在某个物种中，以上一步获取得到的hmm文件为索引，搜寻基因家族成员，以A.hmm.out的结果输出
hmmsearch --domtblout mal.hmm.out B.hmm mal.clean.faa

#(3)blastp和hmm的结果取交集或者取合集，得到每个物种的结果之后，提取相应的氨基酸序列
seqkit grep -f agl agl.clean.faa > agl.PF00089.pep

#(4)先做MSA，此时在tree环境下工作
#在比对建树之前一定要先记得修改序列的名字！！！！！！
~/miniconda3/envs/tree/bin/muscle -align agl.mal.tca.fasta -output agl.mal.tca.align.afa

#(5)用多序列比对结果建树
~/miniconda3/envs/tree/bin/iqtree2 -s agl.mal.tca.align.afa -nt 20



*****GH*****
使用blastx：收集的cDNA序列按照六种方式翻译成蛋白质后检索蛋白质数据库。（不使用基因组作为database是因为结果是scaffold号而不是基因号，后续按照坐标搜索太麻烦，所以直接用蛋白质序列吧）

检索gh1、gh5、gh7、gh28、gh31、gh45、gh48
blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/maldb/mal -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh28.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh28.cdna.blastn.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh1.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh1/agl.gh1.cdna.blastx.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh5.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh5/agl.gh5.cdna.blastx.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh31.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh31/agl.gh31.cdna.blastx.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh28.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh28/agl.gh28.cdna.blastx.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh45.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh45/agl.gh45.cdna.blastx.out

blastx -db /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/agldb/agl -query /malgenome/zfiles/Lepidoptera/9genefamily/gh/gh48.cdna.fasta -num_threads 20 -outfmt 6 -evalue 1e-5 -out /malgenome/zfiles/Lepidoptera/9genefamily/gh/1blastp/out/gh48/agl.gh48.cdna.blastx.out



hmmbuild 1.hmm GH1.seed.ann

hmmsearch --domtblout mal.abc.hmm.out abc.hmm mal.clean.faa

hmmsearch --domtblout mal.gst2.hmm.out gst2.hmm mal.clean.faa
agl agr atu bmo cse dme dpl dva dvi ead lde mal mse ota pma ppy pxy rbi rma tca tni
 




10 RNA-seq
#sratoolkit已经安装在了software里，已经编译到环境变量里了，直接用就行
#下载sra文件
prefetch --option-file SRR_Acc_List.txt
#将sra文件转换为fastq格式
fasterq-dump -O /malgenome/zfiles/Lepidoptera/10RNA-seq/1larva/ -S -p -e 16 SRR26116071.sra --include-technical 
fasterq-dump -O /malgenome/zfiles/Lepidoptera/10RNA-seq/2adult/ -S -p -e 16 SRR26116073.sra --include-technical 
fasterq-dump -O /malgenome/zfiles/Lepidoptera/10RNA-seq/3pupa/ -S -p -e 16 SRR26116072.sra --include-technical
#压缩
gzip SRR26116071_?.fastq 

-----------------------------------------------------------
#压缩完之后先看下数据质量
#使用fastqc和multiqc的时候记得用全路径。它们装在了rnaseq中。
/miniconda3/envs/rnaseq/bin/fastqc -t 20 -o /malgenome/zfiles/Lepidoptera/10RNA-seq/1larva/ SRR26116071_1.fastq.gz > qc.log
/miniconda3/envs/rnaseq/bin/multiqc /malgenome/zfiles/Lepidoptera/10RNA-seq/1larva/*.zip -o /malgenome/zfiles/Lepidoptera/10RNA-seq/1larva/ > multiqc.log
-----------------------------------------------------------
#或者直接查看数据质量和质控一起搞完 fastp也在rnaseq环境里
for i in {789..800}
do
/miniconda3/envs/rnaseq/bin/fastp -i CRR592${i}_1.fastq.gz -o CRR592${i}_1.fq.gz -I CRR592${i}_2.fastq.gz -O CRR592${i}_2.fq.gz -w 20
done

#质控完成之后开始比对，align回松墨天牛基因组上
#samtools和hisat2安装在了ngs环境里
#首先用hisat2构建索引。
/miniconda3/envs/ngs/bin/hisat2-build -p 20 genome.fa genome
#索引构建完成之后开始mapping
for i in 795 796
do
/miniconda3/envs/ngs/bin/hisat2 -p 20 -x /malgenome/zfiles/Lepidoptera/10RNA-seq/5malref/genome -1 CRR592${i}_1.fq.gz -2 CRR592${i}_2.fq.gz -S CRR592${i}.sam
done
#比对完之后转换文件格式
for i in 799 800
do
/miniconda3/envs/ngs/bin/samtools view SRR26116071.sam -b > SRR26116071.bam
/miniconda3/envs/ngs/bin/samtools sort SRR26116071.bam -o SRR26116071_sorted.bam  
/miniconda3/envs/ngs/bin/samtools index SRR26116071_sorted.bam
done
#转换完格式之后先把gff3注释文件转成gft格式，使用gffread，装在rnaseq里
/miniconda3/envs/rnaseq/bin/gffread genome.gff3 -T -o genome.gtf
#开始定量
/miniconda3/envs/rnaseq/bin/featureCounts -T 20 -p -a /malgenome/zfiles/Lepidoptera/10RNA-seq/5malref/genome.gtf -o larva.txt SRR26116071_sorted.bam 





11.正选择基因分析
#此分析中的软件都在paml环境下
#!/bin/bash
#PBS -N mafft
#PBS -l nodes=5:ppn=4
#PBS -q blade
#PBS -l walltime=9999:00:00
#PBS -j oe
#(1)先比对这些单拷贝的OG
cd /malgenome/zfiles/Lepidoptera/11positive/2Single_Copy_Orthologue_Sequences/
for i in {5821..6098}
do
/miniconda3/envs/paml/bin/mafft --thread 20 --maxiterate 1000 --globalpair OG000${i}.fa > OG000${i}.align.fa
done

#(2)从数据库里下载物种的CDS文件
#或者用gffread从genome和gff文件中提取每个物种的CDS文件，在rnaseq环境下，但是不好用，在后续pal2nal的时候报错了
#gffread genes.gff3 -g Pterostichus_madidus-GCA_911728475.2-softmasked.fa -x pma_cds.fasta

#先整理好cds文件，把标题行中多余的信息删除，这里的代码意思是删除掉标题行中出现的第一个空格之后的信息
sed '/^>/s/ .*//' file1 > file2
#将标题行中的氨基酸序列号留下，其他内容删除，这里的其他内容指的是CDS序列的名称。这里代码意思是将>lcl|和_cds_以及之间的内容换成>，[^>]*指的是>以外的任意字符匹配任意次
sed 's/^>lcl|[^>]*_cds_/>/g' file1 > file2
#将序列名称中的最后一个_以及后面的数字删除
sed 's/_[0-9]*$//' file1 > file2
#从数据库中获取的CDS序列文件的序列名称中含有对应的蛋白质序列ID，因此在对CDS序列文件的名称进行整理的时候将其他的信息都删除了，只留下了对应的蛋白质序列的ID。
#有个别物种的CDS序列文件中名称行含有的是CDS序列的ID，不是对应的蛋白质序列的ID，因此这种情况需要利用基因组的gff注释文件来找到CDS和蛋白质序列ID的对应关系。这一步在Excel就可以很方便的进行。

#(3)根据id list和物种的CDS文件提取出相对应的CDS序列。
#注意，在这一步利用的atu.cds.fa CDS文件是已经处理过的、将序列转换成相应氨基酸序列id的文件。而atu.id则是567个单拷贝基因的氨基酸id。
#提取每个物种的567个CDS序列
seqkit grep -f agl.id agl.cds.fa > agl.cds.fasta

#(4)对于rma，要把提取得到的cds序列文件中的cds ID换成相应的氨基酸ID
# 读取rma.cds.aa.id文件并创建一个字典
id_mapping = {}
with open('rma.cds.aa.id', 'r') as id_file:
    for line in id_file:
        cds_id, aa_id = line.strip().split('\t')
        id_mapping[cds_id] = aa_id

# 读取fasta序列文件，进行名称替换并输出
output_lines = []
with open('rma.cds.fasta', 'r') as fasta_file:
    for line in fasta_file:
        if line.startswith('>'):
            cds_id = line[1:].strip()
            aa_id = id_mapping.get(cds_id, cds_id)
            output_lines.append(f'>{aa_id}\n')
        else:
            output_lines.append(line)

# 将结果写入输出文件
with open('rma.aaid.cds.fasta', 'w') as output_file:
    output_file.writelines(output_lines)

#(5)按照pep比对的顺序提取出氨基酸序列的id号，同时也是CDS序列的顺序
cd /malgenome/zfiles/Lepidoptera/11positive/3pepalign/
for i in {5821..7054}
do
awk '/^>/ {print substr($1, 2)}' OG000${i}.aln.fa > OG000${i}.aaid
done
#得到的aaid文件保存到/malgenome/zfiles/Lepidoptera/11positive/4id/中

#(6)按照这个cds的id顺序来提取cds序列
cd /malgenome/zfiles/Lepidoptera/11positive/4id/
#使用bash脚本来完成要求，以下是内容
all_cds_file="all.cds.fa"

# 循环提取序列
for ((i=5821; i<=7054; i++)); do
    idlist_file="OG000${i}.aaid"
    output_file="OG000${i}.cds.fa"
    
    # 清空输出文件
    > "$output_file"

    # 逐行读取idlist文件
    while IFS= read -r id; do
        # 使用awk提取对应id的序列并删除空行
        awk -v RS=">" -v id="$id" '$1 == id {print ">" $0}' "$all_cds_file" | awk 'NF' >> "$output_file"
    done < "$idlist_file"
done
#提取出来的cds序列转移到/malgenome/zfiles/Lepidoptera/11positive/5ogcds/

#将序列cds的序列id转换成相对应物种的名称，使用python脚本完成，同时将pepaln的比对文件的序列id也转换为物种名称。以下是内容
import os
def process_file(input_filename, output_filename, id_mapping):
    with open(input_filename, 'r') as input_file, open(output_filename, 'w') as output_file:
        current_sequence = ''
        for line in input_file:
            if line.startswith('>'):
                current_sequence = line.strip()[1:]
                if current_sequence in id_mapping:
                    output_file.write(f'>{id_mapping[current_sequence]}\n')
                else:
                    output_file.write(line)
            else:
                output_file.write(line)

def main():
    # 读取cds.aa.id文件，创建一个字典用于存储匹配关系
    id_mapping = {}
    with open('aa.species.id', 'r') as id_file:
        for line in id_file:
            cols = line.strip().split('\t')
            if len(cols) == 2:
                id_mapping[cols[0]] = cols[1]

    # 遍历文件并处理每个文件
    input_files = [f'OG000{i}.aln.fa' for i in range(5821, 7055)]  # 修改范围以匹配你的文件名
    for input_filename in input_files:
        # 检查文件是否存在
        if os.path.exists(input_filename):
            output_filename = f'{input_filename[:-6]}.aln.speciesid.fa'  # 生成新的输出文件名
            process_file(input_filename, output_filename, id_mapping)
        else:
            print(f"Warning: File {input_filename} not found. Skipping.")

if __name__ == "__main__":
    main()

#(7)#使用pal2nal将转换好的比对文件转换成密码子的比对
#这个命令可以将path1下的空文件全部转移到path2中，记住要全路径
find path1 -type f -empty -exec mv {} path2 \;

#开始pal2nal
cd /malgenome/zfiles/Lepidoptera/11positive/6pal2nal/
for i in {5821..7054}
do
perl /miniconda3/envs/paml/bin/pal2nal.pl OG000${i}.oneline.pepaln.fa OG000${i}.oneline.cds.fa -output paml -nogap > OG000${i}.codon
done

#(8)开始codeml
#####null model
cd /malgenome/zfiles/Lepidoptera/11positive/6codeml0/

start_index=5821
end_index=6282

# 循环处理每个文件
for ((i=start_index; i<=end_index; i++)); do
    # 构建文件名
    file_name="OG000${i}.codon"
    
    # 替换 codeml.ctl 中的相关内容
    sed -i "s/ALN/${file_name}/" codeml01.ctl
    sed -i "s/TREE/species.tree/" codeml01.ctl
    sed -i "s/OUT/${file_name%.codon}.null.txt/" codeml01.ctl  # 去除文件扩展名 ".codon"
    
    # 运行 codeml
    /miniconda3/envs/paml/bin/codeml codeml01.ctl
    
    # 恢复 codeml.ctl 的原始状态，以备下一次循环使用
    cp codeml.0.ctl.bck codeml01.ctl
done

#####codeml.0.ctl.bck的内容
      seqfile = ALN            * Path to the alignment file
     treefile = TREE           * Path to the tree file
      outfile = OUT            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = 1           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = 2         * Models for ω varying across lineages
	  NSsites = 2          * Models for ω varying across sites
    CodonFreq = 7        * Codon frequencies
	  estFreq = 0        * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
    fix_omega = 1         * Estimate or fix omega
        omega = 1        * Initial or fixed omega


#####alternative model
cd /malgenome/zfiles/Lepidoptera/11positive/6codemlal/

start_index=5821
end_index=6282

# 循环处理每个文件
for ((i=start_index; i<=end_index; i++)); do
    # 构建文件名
    file_name="OG000${i}.codon"
    
    # 替换 codeml.ctl 中的相关内容
    sed -i "s/ALN/${file_name}/" codemlal1.ctl
    sed -i "s/TREE/species.tree/" codemlal1.ctl
    sed -i "s/OUT/${file_name%.codon}.al.txt/" codemlal1.ctl  # 去除文件扩展名 ".codon"
    
    # 运行 codeml
    /miniconda3/envs/paml/bin/codeml codemlal1.ctl
    
    # 恢复 codeml.ctl 的原始状态，以备下一次循环使用
    cp codeml.al.ctl.bck codemlal1.ctl
done

#####codeml.al.ctl.bck的内容
      seqfile = ALN            * Path to the alignment file
     treefile = TREE           * Path to the tree file
      outfile = OUT            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = 1           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = 2         * Models for ω varying across lineages
	  NSsites = 2          * Models for ω varying across sites
    CodonFreq = 7        * Codon frequencies
	  estFreq = 0        * Use observed freqs or estimate freqs by ML
        clock = 0          * Clock model
    fix_omega = 0         * Estimate or fix omega
        omega = 1        * Initial or fixed omega

#(9)开始统计
#将null和alternative的codeml结果中的lnL提取到同一个文件中，第一列是文件名，第二列是提取出来的lnL(似然值的对数值)
#null model
# 循环处理文件
for i in {5821..7054}; do
    filename="OG000${i}.null.txt"
    
    # 提取文件名并写入结果文件
    echo -n "${filename}" >> lnL.branchsite.null
    
    # 提取lnL内容并写入结果文件
    grep 'lnL' "${filename}" | sed 's/..*\:\ *//' | sed 's/\ ..*//' | tr '\n' '\t' >> lnL.branchsite.null
    
    # 添加换行符
    echo "" >> lnL.branchsite.null
done

#alternative model
# 循环处理文件
for i in {5821..7054}; do
    filename="OG000${i}.al.txt"
    
    # 提取文件名并写入结果文件
    echo -n "${filename}" >> lnL.branchsite.al
    
    # 提取lnL内容并写入结果文件
    grep 'lnL' "${filename}" | sed 's/..*\:\ *//' | sed 's/\ ..*//' | tr '\n' '\t' >> lnL.branchsite.al
    
    # 添加换行符
    echo "" >> lnL.branchsite.al
done

#(10)计算alternative的lnL与null的lnL之间差值的两倍，公式是：似然比（Likelihood Ratio，LR），即两个模型似然度的比值。LR = 2 * (ln(likelihood(alternative model)) - ln(likelihood(null model)))。
#这里在Excel做就可以，直接相减，因为从OG000${i}.al.txt和OG000${i}.null.txt中提取出来的那个是对数后的值。
完成之后进入R
