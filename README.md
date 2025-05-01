# SAPTyper_v1.0


## SAPTyper软件安装和使用
1. 主体软件不需要额外安装，只需要直接将整个目录拷贝到某个文件夹下即可。建议在D盘根目录建立SAPTyper的目录，将整个文件夹下的内容直接拷贝过去。

2. 需要预安装如下软件：
a. 需要预先安装perl，dependence目录下提供了一个64位的perl软件可直接安装
b. 需要预先安装pFind或MaxQuant。根据项目测试，建议安装pFind。pFind需要申请license，请访问“http://pfind.ict.ac.cn/”获得最新的pFind软件

## 3. 配置：

a. 需要配置bin\parameter目录下的以下文件：
mqpar.xml：Maxquant的参数文件，用Maxquant设置好参数后测试搜索，即可得到这文件
pFind.cfg：pFind的参数文件，用pFind设置好参数后测试搜索，即可得到这文件
pParse.cfg：pFind的参数文件，用pFind设置好参数后测试搜索，即可得到这文件
default.parameter.par：需要设置public_var、ref_parameter、pparse_parameter、seq_db、SNP_Anno、SAP_Anno的路径（根据实际安装路径，在bin的data目录下）

b. 需要配置bin目录下的文件：
dep_program.dir：设置Maxquant和pFind的安装路径，到bin目录；
anno_data_file.txt：注释文件的路径，一般在bin下的data里

## 4. 运行：

a. 运行参数设置：
把bin下的default目录下的三个文件，拷贝到运行目录中（数据可以不在这目录下，但需要设置参数文件指向数据）
根据实际数据，设置这三个文件相关参数

b. 在运行目录下，执行命令行："绝对路径\SAP.analyse.pl project.par"即可
