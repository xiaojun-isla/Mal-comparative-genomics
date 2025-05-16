/Lepidoptera/1genome/

这个文件夹下主要是在提取蛋白质序列的ID。

原材料文件：
genomic.gff3  基因组注释文件
pep.fa        物种的蛋白质序列文件，在TBtools中用的
id.txt        利用TBtools提取到的代表性转录本的ID，也就是CDS的ID

产生的文件：
output.txt         从gff文件中提取出来的第三列为CDS的行的第9列
output_unique.txt  output.txt文件排序、去重后的文件
final_id.txt       联合output_unique.txt、id.txt两个文件提取出来的蛋白质ID
final_ids.txt      排序、去重后的最后的蛋白质序列ID，用这个文件中的ID去TBtools中提取蛋白质序列

使用TBtools提取出来的ID是CDS序列的ID，并不是蛋白质序列的ID，这些CDS序列的ID来源于gff文件。
因此若想得到CDS对应的蛋白质序列，需要先得到CDS对应的蛋白质序列ID，而蛋白质序列的ID从何而来呢？
答：从gff文件的第9列而来。

首先使用TBtools来从gff文件中获得代表性转录本的ID，也就是CDS的ID，oriid.txt整理之后就是用于从gff文件中查找蛋白质序列ID的CDS的id.txt文件。
从gff文件中提取出第三列为CDS的行的第九列，得到文件output.txt，然后排序、去重，得到output_unique.txt，用id.txt中的转录本ID来搜查output_unique.txt中的匹配行，然后提取出protein_id部分，也有可能是ID=cds部分，视具体情况而定，得到final_id.txt，然后排去重，得到最终的protein_ID，即final_ids.txt文件。





/Lepidoptera/2ortho/
将上一步得到的蛋白质序列文件整理好，用统一格式命名，以faa结尾，然后放在一个文件夹里，使用orthofinder时制定这个为输入文件的文件夹。

蛋白质序列文件也可以以.fa .faa .fasta .fas .pep结尾





/Lepidoptera/3cafe/
原材料文件：
SpeciesTreeAlignment.fa
SpeciesTree_rooted.txt
Orthogroups.GeneCount.tsv

生成文件：
准备cafe需要的输入文件之一，一个树文件，带有分化时间：
r8s_ctl_file.txt  
r8s_tmp.txt
r8s_ultrametric.txt 准备好的超度量树
准备cafe需要的输入文件之一，一个基因家族数目统计文件：
cafe.input.tsv
filtered.cafe.input.tsv  准备好的基因家族统计文件





/Lepidoptera/4phylogeny/
系统发育分析需要的原始软件：
SpeciesTreeAlignment.fa  来源于OrthoFinder/Results_Jul10/MultipleSequenceAlignments/，根据教程，该文件为orthofinder生成好的单拷贝基因的比对文件

产生的中间文件及结果文件的注释：
SpeciesTreeAlignment_trimal.fa
RAxML_bestTree.out
r8s_ctl_file.txt
r8s_tmp.txt
r8s_ultrametric.txt 最终产生的超度量树

输出文件解读：
在这里选用了赤拟谷盗作为外群
RAxML_bestTree.out                    最终选择的最佳进化树
RAxML_bipartitions.out                加入了自展值的进化树，自展值标在小括号外面
RAxML_bipartitionsBranchLabels.out    与RAxML_bipartitions.out一样，只不过自展值标在方括号内，一般选择该文件进行后续进化树美化，要在mega中打开的话，要将out后缀改为nwk
RAxML_bootstrap.out                   每次bootstrap分析建出来的进化树结构，共100次分析，故文件中包含100个进化树
RAxML_info.out                        日志信息





/Lepidoptera/5GOKEGG/
使用到的原始文件：
Gamma_change.tab
Gamma_count.tab
Gamma_family_results.txt
Orthogroups.tsv
Orthogroups.txt
all.faa

产生的文件：
p0.05.significant.ogs 0.05显著收缩或扩张的OG ID
p0.01.significant.ogs 0.01显著收缩或扩张的OG ID
Gamma_p0.05change.tab 显著扩张 / 收缩的基因家族在每个节点的收缩与扩增数目
Gamma_p0.01change.tab 显著扩张 / 收缩的基因家族在每个节点的收缩与扩增数目
mal.expanded 松墨节点上发生显著扩张的基因家族，0.05
mal.contracted 松墨节点上发生显著收缩的基因家族，0.05
mal.0.05.expanded.genes 松墨0.05显著扩张的gene ID
mal.0.05.contracted.genes 松墨0.05显著收缩的gene ID
mal.0.01.expanded.genes 松墨0.01显著扩张的gene ID
mal.0.01.contracted.genes 松墨0.01显著收缩的gene ID
mal.0.05.expanded.pep.fa 0.05显著扩张的蛋白序列
mal.0.01.expanded.pep.fa 0.01显著扩张的蛋白序列
mal.0.05.contracted.pep.fa 0.05显著收缩的蛋白序列
mal.ogs 含有松墨基因的OG ID(背景基因的OG ID)
mal.ogs.genes 含有松墨基因的OG ID所包含的所有松墨基因
mal.ogs.pep.fa 含有松墨基因的OG ID所包含的所有松墨基因的相应的序列(背景基因的序列)
go-basic.obo 原始注释库
GO.library 修饰好的注释库
GO.id 注释库的ID号
GO.name 注释库中的名字
GO.class 注释库的三大本体论词条


/Lepidoptera/8synteny/
MCScanX的结果不太理想，结果比较杂乱，不方便观察，也不能显示出最想要显示的部分结果，因此使用MCScan

原材料文件：
gff文件
cds文件
