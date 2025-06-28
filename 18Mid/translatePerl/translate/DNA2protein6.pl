#!/usr/bin/perl
use warnings;
use strict;




###################################################################################
#此脚本与DNA2protein5.pl的差别在于用来做blast的库分别做了在96%水平的去冗余

#输入流感病毒的DNA序列，输出为翻译好的蛋白序列
#蛋白序列在基因序列名字后面加上蛋白名字，特别是对于MP和NS蛋白
#基本思路如下：
#首先，对于每条DNA序列做blastx,得到它的蛋白序列和所属基因；
#其次，如果某条DNA序列属于编码两条蛋白的序列，如MP,NS,PA,PB1-F2,则
#继续用该DNA序列对于另外一个只包含M2，NS2，PA-X和PB1-F2的库做blastX，得到
#它的编码序列
#输入参数如下：
#$infile   需要翻译的DNA序列文件名（完整路径）
#$blastDir 安装blast工具的目录，需要在最后加“/"
#$blastLibraryDir  #压缩包中的forblast文件夹
#$tempDir   #临时文件夹目录
#$referSeqDir   ##压缩包中的data文件夹
#$proteinSeqOutfile  #输出文件

#运行示例如下：
#perl DNA2protein6.pl ./temp/DNAseq /usr/bin/ ./forblast/ ./temp/ ./data/ ./temp/proteinSeq
####################################################################################
my ($infile,$blastDir,$blastLibraryDir,$tempDir,$referSeqDir,$proteinSeqOutfile,$mafftDir)=@ARGV;
if(scalar @ARGV < 6){
	die"error!Please input the DNA sequence file,blastDir,blastLibraryDir,tempDir,referSeqDir and the outfile for protein sequences\n";
}
my $start=time();


#首先，读取需要翻译的DNA序列
open(IN,$infile)or die"$!";
my @allid;
my $seqRef;
my $id;
my $new2old;
my $newSeqNum=1;
while(<IN>){
	s/^[\r\t\s\n]{1,}//;
	s/[\r\t\s\n]{1,}$//;
	chomp;
	next if(/^$/);
	if(/^>(.+)$/){
		my $oldID=$1;
		$id="Userseq$newSeqNum";
		$new2old->{$id}=$oldID;
		$newSeqNum++;
		push(@allid,$id);
	}else{
		s/-/N/g;  #处理序列中包含gap的情况
		my @line=split(//,uc($_));
		push(@{$seqRef->{$id}},@line);
	}
}
close IN or die"$!";


#其次，对于不包含M2，NS2，PB1-F2和PA-X的库做blastx
my $DNAseqFile=$tempDir."DNAseq";
open(OUT,">$DNAseqFile")or die"$!";
foreach my$id(keys %{$seqRef}){
	print OUT ">$id\n",join("",@{$seqRef->{$id}}),"\n";
}
close OUT or die"$!";

my $database1=$blastLibraryDir."blastLibrary1";
my $blastOutfile1=$tempDir."blastOut1";

`$blastDir"blastall" -p blastx -i $DNAseqFile -d $database1 -o $blastOutfile1 -a 2 -F F -v 1 -b 1`;
#`"blastx" -query $DNAseqFile -db $database1 -out $blastOutfile1 -num_threads 2 -seg no -max_target_seqs 1`;

my ($id2seq1,$id2info1)=getProteinFromBlastx($blastOutfile1);
#id2seq->{$id}=$oneseq
#id2info->{$id}->{'target'}=$target (ID_gene_type)
#id2info->{$id}->{'start'}=$start
#id2info->{$id}->{'end'}=$end

#接下来，对于输入序列中的M，NS，PB1和PA基因，继续对第二个库做blastx
my @remain1;
my @remain2;
foreach my$id(keys %{$id2info1}){
	my $target=$id2info1->{$id}->{'target'};
	my @targetName=split(/\_/,$target);
	my $gene=$targetName[-2];
	my $type=$targetName[-1];

	#如果是A型流感的M，NS和PA基因，或者B/C型的NS基因，则通过多序列比对的方法来确定第二个编码区
	if( ($type eq 'A' && ($gene=~/M/i or $gene=~/NS/i or $gene=~/PA/i) ) or ($type =~/[BC]/i && $gene=~/NS/i) ){
		push(@remain1,$id);

	#如果是A型流感的PB1基因，或者B/C型的M基因，则继续通过blast的方法确定第二个编码区
	}elsif( ($type =~/[BC]/i && $gene=~/M/i) or ($type eq 'A' && $gene=~/PB1/i) ){
		push(@remain2,$id);
	}
}

#如果存在A型流感的PB1基因，或者B/C型的M基因，则继续通过blast的方法确定第二个编码区
my ($id2seq2,$id2info2);
if(scalar @remain2 > 0){
	my $DNAseqFile_remain2=$tempDir."DNAseq_remain2";
	open(OUT,">$DNAseqFile_remain2")or die"$!";
	foreach my$id(@remain2){
		print OUT ">$id\n",join("",@{$seqRef->{$id}}),"\n";
	}
	close OUT or die"$!";
	
	my $database2=$blastLibraryDir."blastLibrary2";
	my $blastOutfile2=$tempDir."blastOut2";

#    `"blastx" -query $DNAseqFile_remain2 -db $database2 -out $blastOutfile2 -num_threads 2 -seg no -max_target_seqs 1`;

	`$blastDir"blastall" -p blastx -i $DNAseqFile_remain2 -d $database2 -o $blastOutfile2 -a 2 -F F -v 1 -b 1`;

	($id2seq2,$id2info2)=getProteinFromBlastx($blastOutfile2);
}

#如果存在A型流感的M，NS和PA基因，或者B/C型的NS基因，则通过多序列比对的方法来确定第二个编码区，然后直接翻译出来
my $id2seq3;
my $id2info3;
if(scalar @remain1 > 0){
	my $geneRef; #先对基因进行分类
	foreach my$id(@remain1){
		my $target=$id2info1->{$id}->{'target'};
		my @targetName=split(/\_/,$target);
		my $gene=$targetName[-2];
		my $type=$targetName[-1];
		push(@{$geneRef->{$type}->{$gene}},$id);
	}
	foreach my$type(sort keys %{$geneRef}){
		foreach my$gene(sort keys %{$geneRef->{$type}}){
			my @id=@{$geneRef->{$type}->{$gene}};
			my $DNAseqFile_remain1=$tempDir."DNAseq_remain1";
			open(OUT,">$DNAseqFile_remain1")or die"$!";
			foreach my$id(@id){
				print OUT ">$id\n",join("",@{$seqRef->{$id}}),"\n";
			}
			close OUT or die"$!";

			my $DNAseqFile_remain1_align=$tempDir."DNAseq_remain1_align";
			my $geneNew=$gene; #为了处理NS1，NS2，M1，M2，PA-X的问题，它们对应的基因分别是NS，M和PA，因此需要处理后面的后缀，方便读取参考序列文件
			$geneNew=~s/[0-9]{1,}$//ig;
			$geneNew=~s/\-X$//ig;
			my $aligned2original=align2referSeq2($DNAseqFile_remain1,$type,$geneNew,$referSeqDir,$tempDir,$DNAseqFile_remain1_align,$mafftDir); #把用户的序列变成与参考序列一样
			my ($secondCodingRef,$codingRegionRef)=getSecondCodingRegion2($DNAseqFile_remain1_align,$type,$gene,$aligned2original);#获取第二个编码区序列
			foreach my$id(keys %{$secondCodingRef}){
				my $codingSeq=join("",@{$secondCodingRef->{$id}});
				my $proteinSeq=DNA2protein($codingSeq,0);
				$proteinSeq=~s/^\-{1,}//;  #去掉前面和后面的gap
				$proteinSeq=~s/\-{1,}$//;
				$id2seq3->{$id}=$proteinSeq;
				$id2info3->{$id}=$codingRegionRef->{$id};
			}
		}
	}
}

#最后，输出所有蛋白序列，在每个序列名字后面加上其编码的蛋白
open(OUT,">$proteinSeqOutfile")or die"$!";
foreach my$id(keys %{$id2seq1}){
	my $target=$id2info1->{$id}->{'target'};
	my @targetName=split(/\_/,$target);
	my $gene=$targetName[-2];
	my $type=$targetName[-1];
	
	my $old=$new2old->{$id};
	my $geneName=$gene;
	if($gene=~/M/i){
		if($type=~/A/i){
			$geneName="M1";
		}elsif($type=~/B/i){
			$geneName="BM1";
		}elsif($type=~/C/i){
			$geneName="CM1";
		}
	}elsif($gene=~/NS/i){
		$geneName="NS1";
	}
	my $seqName=$old."_".$geneName;
	my $codingStart1=$id2info1->{$id}->{'start'};  #增加编码区对应的区域
	my $codingEnd1=$id2info1->{$id}->{'end'};
	$seqName.="(".$codingStart1."..".$codingEnd1.")";

	unless($id2seq1->{$id} eq ''){    #需要考虑序列为空的情况
		print OUT ">$seqName\n",$id2seq1->{$id},"\n";
	}
	if(exists $id2seq2->{$id}){
		my $geneName=$gene;
		if($gene=~/M/i){
			if($type=~/B/i){
				$geneName="BM2";
			}elsif($type=~/C/i){
				$geneName="CM2";
			}
		}elsif($gene=~/PB1/i){
			$geneName="PB1-F2";
		}
		my $seqName=$old."_".$geneName;

		my $codingStart2=$id2info2->{$id}->{'start'};  #增加编码区对应的区域
		my $codingEnd2=$id2info2->{$id}->{'end'};
		$seqName.="(".$codingStart2."..".$codingEnd2.")";

		unless($id2seq2->{$id} eq ''){         #需要考虑序列为空的情况
			print OUT ">$seqName\n",$id2seq2->{$id},"\n";
		}
	}

	if(exists $id2seq3->{$id}){
		my $geneName=$gene;
		if($gene=~/M/i){
			$geneName="M2";
		}elsif($gene=~/PA/i){
			$geneName="PA-X";
		}elsif($gene=~/NS/i){
			$geneName="NS2";
		}
		my $seqName=$old."_".$geneName;

		my $codingRegion=$id2info3->{$id};  #增加编码区对应的区域
		$seqName.="(".$codingRegion.")";

		unless($id2seq3->{$id} eq ''){    #需要考虑序列为空的情况
			print OUT ">$seqName\n",$id2seq3->{$id},"\n";
		}
	}
}
close OUT or die"$!";

my $end=time();
my $timeConsume=$end-$start;
print "Job done. Time consumes: $timeConsume\n";








############################################################
#此函数与getSecondCodingRegion的区别在于给出了第二个编码区的位置

#此函数是对于已经比对到参考序列的序列，获取其第二个编码区的序列
#该函数只适合A型的M，NS和PA，以及B和C的NS基因
#输入为比对到参考序列的序列文件，以及流感型和基因
#输出为每条序列第二个编码区的序列以及对应的原始序列的位置：$codingRef->{$id}=\@oneseq
#$codingRegionRef->{$id}=XX..XX
###########################################################
sub getSecondCodingRegion2{
	my ($DNAseqFile_remain1_align,$type,$gene,$aligned2original)=@_;
	
	my $seqRef=readSeq($DNAseqFile_remain1_align);
	my $codingRef;
	my $codingRegionRef;
	foreach my$id(keys %{$seqRef}){
		my @oneseq=@{$seqRef->{$id}};
		my @coding;
		my @pos;
		if($type =~/A/i){
			if($gene=~/M/i){
				push(@coding,@oneseq[0..25,714..981]);
				my $secondExon=join("",@oneseq[714..981]);  #如果序列只包含第一个外显子，则只考虑第一个外显子对应的位置
				for(my $i=0;$i<=981;$i++){
					next if($i>25 && $i<714);
					last if($i > 25 && ($secondExon!~/[atcg]/i));  #如果序列只包含第一个外显子，则只考虑第一个外显子对应的位置
					my $pos=$i+1;
					my $originalPos=$aligned2original->{$id}->{$pos};
					push(@pos,$originalPos);
				}
			}elsif($gene=~/NS/i){
				push(@coding,@oneseq[0..29,502..837]);
				my $secondExon=join("",@oneseq[502..837]);
				for(my $i=0;$i<=837;$i++){
					next if($i>29 && $i<502);
					last if($i > 29 && ($secondExon!~ /[atcg]/i));
					my $pos=$i+1;
					my $originalPos=$aligned2original->{$id}->{$pos};
					push(@pos,$originalPos);
				}
			}elsif($gene=~/PA/i){
				push(@coding,@oneseq[0..569,571..759]);
				my $secondExon=join("",@oneseq[571..759]);
				for(my $i=0;$i<=759;$i++){
					next if($i>569 && $i<571);
					last if($i > 569 && ($secondExon!~ /[atcg]/i));
					my $pos=$i+1;
					my $originalPos=$aligned2original->{$id}->{$pos};
					push(@pos,$originalPos);
				}
			}
		}elsif($type=~/B/i){
			if($gene=~/NS/i){
				push(@coding,@oneseq[0..35,691..1026]);
				my $secondExon=join("",@oneseq[691..1026]);
				for(my $i=0;$i<=1026;$i++){
					next if($i>35 && $i<691);
					last if($i > 35 && ($secondExon!~ /[atcg]/i));
					my $pos=$i+1;
					my $originalPos=$aligned2original->{$id}->{$pos};
					push(@pos,$originalPos);
				}
			}
		}elsif($type=~/C/i){
			if($gene=~/NS/i){
				push(@coding,@oneseq[0..186,500..861]);
				my $secondExon=join("",@oneseq[500..861]);
				for(my $i=0;$i<=861;$i++){
					next if($i>186 && $i<500);
					last if($i > 186 && ($secondExon!~ /[atcg]/i));
					my $pos=$i+1;
					my $originalPos=$aligned2original->{$id}->{$pos};
					push(@pos,$originalPos);
				}
			}
		}
		push(@{$codingRef->{$id}},@coding);

		my $region=number2region(@pos);
		$codingRegionRef->{$id}=$region;
	}

	return ($codingRef,$codingRegionRef);
}


############################################################################################
#此程序是把一串数字（单调非递减）转化为基因编码区域。比如1,2,3,5,6,8,9转化为1..3,5..6,8..9
#同时需要考虑密码子的编码情况，比如起点和终点只能位于能够整除3的位置
############################################################################################
sub number2region{
	my @number=@_;

	my ($start,$end)=(0,0);
	my @region;
	for(my $i=0;$i<scalar @number;$i++){
		next if($number[$i]<=0);
		if((scalar @region)==0 && ($i+1)%3 ==1){
			$start=$number[$i];
			push(@region,$start);
		}
		if($i >0){
			next if((scalar @region)==0);
			if( ($number[$i] - $number[$i-1]) > 1){
				$end=$number[$i-1];
				push(@region,$end);
				$start=$number[$i];
				push(@region,$start);
			}elsif($i==(scalar @number-1)){
				my $k=$i;
				while(($number[$k]-$number[$k-1]) == 0){
					$k--;
				}
				if(($k+1)%3 ==0){
					$end=$number[$k];
				}else{
					if(($k+1)%3 ==1){
						$end=$number[$k-1];
					}elsif(($k+1)%3==2){
						$end=$number[$k-2];
					}
				}
				push(@region,$end);
			}
		}
	}

	my @pair;
	for(my $i=0;$i<scalar @region;$i+=2){
		next if($region[$i] > $region[$i+1]);
		push(@pair,$region[$i]."..".$region[$i+1]);
	}

	return join(",",@pair);
}




##################################################################################################
#此函数是把用户的序列整理成与参考序列相似的形式
#输入为需要处理的序列文件，流感型，基因，参考序列所在文件夹，临时文件夹，以及输出文件名
#输出为标准序列与原有序列每个位点之间的对应关系
##################################################################################################
use Data::Dumper;

sub align2referSeq2{
	my ($DNAseqFile_remain1,$type,$gene,$libraryDir,$tempDir,$DNAseqFile_remain1_align,$mafftDir)=@_;
	
	my $referSeqFile=$libraryDir.$type."_".$gene.".fas";
	my $combinedFile=$tempDir."combineSeq";
	my $combinedFile_align=$combinedFile."_align";
	`cat $DNAseqFile_remain1 $referSeqFile > $combinedFile`;
	`mafft --retree 1 $combinedFile 2>/dev/null > $combinedFile_align`;
#	`$mafftDir"mafft" --retree 1 $combinedFile > $combinedFile_align`;

	my $combinedSeqRef=readSeq($combinedFile_align);
	my $referSeqRef=readSeq($referSeqFile);
	my @referID=keys %{$referSeqRef};
#    my $combinedSeqRef = readSeq($combinedFile_align);
#    print "Combined Sequences:\n", Dumper($combinedSeqRef);
#
#    my $referSeqRef = readSeq($referSeqFile);
#    print "Reference Sequences:\n", Dumper($referSeqRef);
#
#    my @referID = keys %{$referSeqRef};
#    print "Reference IDs:\n", Dumper(\@referID);


	my $siteRef_standard;
	foreach my$id(@referID){
	    $id =~ s/\s+$//;  # 去除行尾的空白字符，包括换行符 #就是这个导致我一直报错！！！Can't use an undefined value as an ARRAY reference at /usr/local/lib/python3.6/site-packages/18Mid/translatePerl/translate/DNA2protein6.pl line 406
		my @oneseq=@{$combinedSeqRef->{$id}};
		for(my $i=0;$i<scalar @oneseq;$i++){
			$siteRef_standard->{$i+1}->{$oneseq[$i]}++;
		}
	}

	my @remainPos;
    foreach my$site(sort {$a<=>$b}keys %{$siteRef_standard}){
		my @key=sort {$siteRef_standard->{$site}->{$b} <=> $siteRef_standard->{$site}->{$a}}keys %{$siteRef_standard->{$site}};
        if($key[0] eq '-'){  #如果标准序列的某个位点是空格，那么该位置应该去掉，除非一些特殊的位点
			if((scalar @remainPos)>=6 && (scalar @remainPos) <=8){
				if($type =~/B/i && $gene =~/NS/i){
					push(@remainPos,$site);
				}
			}
			next;
		}else{
			push(@remainPos,$site);
		}
    }

	#finally, get the aligned sample sequences
	my $aligned2original; #给出处理之后的标准序列中每个位置与原有位置的对应关系
	open(OUT,">$DNAseqFile_remain1_align")or die"$!";
	foreach my$id(keys %{$combinedSeqRef}){
		next if(exists $referSeqRef->{$id});    #不考虑标准序列了；
		my @oneseq=@{$combinedSeqRef->{$id}};
		my @finalSeq;
		my $originalPos=0;  #记录原始序列中每个位置的编号
		for(my $i=0;$i<scalar @oneseq;$i++){
			my $pos=$i+1;
			if($oneseq[$i] ne '-'){
				$originalPos++;
			}

			if(grep $_ eq $pos,@remainPos){
				push(@finalSeq,$oneseq[$i]);
				$aligned2original->{$id}->{(scalar @finalSeq)}=$originalPos;  #得到比对之后得到的标准序列中每个位置与原始序列每个位置的关系
			}
		}
		print OUT ">$id\n",join("",@finalSeq),"\n";
	}
	close OUT or die"$!";
	
	return $aligned2original;
}

############################################################
#此函数是对于已经比对到参考序列的序列，获取其第二个编码区的序列
#该函数只适合A型的M，NS和PA，以及B和C的NS基因
#输入为比对到参考序列的序列文件，以及流感型和基因
#输出为每条序列第二个编码区的序列：$codingRef->{$id}=\@oneseq
###########################################################
sub getSecondCodingRegion{
	my ($DNAseqFile_remain1_align,$type,$gene)=@_;
	
	my $seqRef=readSeq($DNAseqFile_remain1_align);
	my $codingRef;
	foreach my$id(keys %{$seqRef}){
		my @oneseq=@{$seqRef->{$id}};
		my @coding;
		if($type =~/A/i){
			if($gene=~/M/i){
				push(@coding,@oneseq[0..25,714..981]);
			}elsif($gene=~/NS/i){
				push(@coding,@oneseq[0..29,502..837]);
			}elsif($gene=~/PA/i){
				push(@coding,@oneseq[0..569,571..759]);
			}
		}elsif($type=~/B/i){
			if($gene=~/NS/i){
				push(@coding,@oneseq[0..35,691..1026]);
			}
		}elsif($type=~/C/i){
			if($gene=~/NS/i){
				push(@coding,@oneseq[0..186,500..861]);
			}
		}
		push(@{$codingRef->{$id}},@coding);
	}

	return $codingRef;
}


##################################################################################################
#此函数是把用户的序列整理成与参考序列相似的形式
#输入为需要处理的序列文件，流感型，基因，参考序列所在文件夹，临时文件夹，以及输出文件名
#没有输出
##################################################################################################

sub align2referSeq{
	my ($DNAseqFile_remain1,$type,$gene,$libraryDir,$tempDir,$DNAseqFile_remain1_align,$mafftDir)=@_;
	
	my $referSeqFile=$libraryDir.$type."_".$gene.".fas";
	my $combinedFile=$tempDir."combineSeq";
	my $combinedFile_align=$combinedFile."_align";
	`cat $DNAseqFile_remain1 $referSeqFile > $combinedFile`;
#	`$mafftDir"mafft" --retree 1 $combinedFile > $combinedFile_align`;
	`mafft --retree 1 $combinedFile 2>/dev/null > $combinedFile_align`;

	my $combinedSeqRef=readSeq($combinedFile_align);
	my $referSeqRef=readSeq($referSeqFile);
	my @referID=keys %{$referSeqRef};

	my $siteRef_standard;
	foreach my$id(@referID){
		my @oneseq=@{$combinedSeqRef->{$id}};
		for(my $i=0;$i<scalar @oneseq;$i++){
			$siteRef_standard->{$i+1}->{$oneseq[$i]}++;
		}
	}

	my @remainPos;
    foreach my$site(sort {$a<=>$b}keys %{$siteRef_standard}){
		my @key=sort {$siteRef_standard->{$site}->{$b} <=> $siteRef_standard->{$site}->{$a}}keys %{$siteRef_standard->{$site}};
        if($key[0] eq '-'){  #如果标准序列的某个位点是空格，那么该位置应该去掉，除非一些特殊的位点
			if((scalar @remainPos)>=6 && (scalar @remainPos) <=8){
				if($type =~/B/i && $gene =~/NS/i){
					push(@remainPos,$site);
				}
			}
			next;
		}else{
			push(@remainPos,$site);
		}
    }

	#finally, get the aligned sample sequences
	open(OUT,">$DNAseqFile_remain1_align")or die"$!";
	foreach my$id(keys %{$combinedSeqRef}){
		next if(exists $referSeqRef->{$id});    #不考虑标准序列了；
		my @oneseq=@{$combinedSeqRef->{$id}};
		my @finalSeq;
		for(my $i=0;$i<scalar @oneseq;$i++){
			my $pos=$i+1;
			if(grep $_ eq $pos,@remainPos){
				push(@finalSeq,$oneseq[$i]);
			}
		}
		print OUT ">$id\n",join("",@finalSeq),"\n";
	}
	close OUT or die"$!";
}


###################################################
#从blast的结果中抽提出蛋白序列以及编码区和hit
#输入为blast的结果文件
#输出为每条查询DNA序列的蛋白序列，以及蛋白target和基因编码区位置
#id2seq->{$id}=$oneseq
#id2info->{$id}->{'target'}=$target
#id2info->{$id}->{'start'}=$start
#id2info->{$id}->{'end'}=$end
###################################################
sub getProteinFromBlastx{
	my ($blastOutFile)=@_;

	my $id2seq;
	my $id2info;
	my $oneseq="";
	my $tag=0;  #表示只读取同一个hit的一个序列
	open(IN,$blastOutFile)or die"$!";
	while(<IN>){
		chomp;
		if(/^Query\=(.+)$/){
			my $newid=$1;
			if(length($oneseq) > 10){
				$id2seq->{$id}=$oneseq;
			}
			$id=$newid;
			$id=~s/^\s{1,}//;
			$id=~s/\s{1,}$//;
			$oneseq="";
			$tag=0;
		}elsif(/^>(.+)$/ && $tag==0){
			my $target=$1;
			$target=~s/^\s{1,}//;
			$target=~s/\s{1,}$//;
			$id2info->{$id}->{'target'}=$target;
			$tag=1;
		}elsif(/^Query\:/ && $tag==1){
			my @line=split(/\s{1,}/);
			unless(exists $id2info->{$id}->{'start'}){  #如果存在一个query对应的target中有两个区域都匹配，则需要考虑这个问题
				$id2info->{$id}->{'start'}=$line[1];
			}
			$id2info->{$id}->{'end'}=$line[3];
			$oneseq.=$line[2];
			while(<IN>){
				chomp;
				my $first=$_;
				my $second=<IN>;
				chomp $second;
				if(($first=~/^$/ && $second=~/^$/) or $first=~/^\s{1,}score/i or $second=~/^\s{1,}score/i){
					$tag=0;last;
				}
				my $content="NA";
				if($first=~/^Query\:/){
					$content=$first;
				}elsif($second=~/^Query\:/){
					$content=$second;
				}
				if($content=~/^query\:/i){
					my @line=split(/\s{1,}/,$content);
					$id2info->{$id}->{'end'}=$line[3];
					$oneseq.=$line[2];
				}
			}
		}
	}
	close IN or die"$!";

	if(length($oneseq) > 10){#处理最后一条记录
		$id2seq->{$id}=$oneseq;
	}

	return ($id2seq,$id2info);
}



###################################################
#从blast的结果中抽提出蛋白序列以及编码区和hit
#输入为blast的结果文件
#输出为每条查询DNA序列的蛋白序列，以及蛋白target和基因编码区位置
#id2seq->{$id}=$oneseq
#id2info->{$id}->{'target'}=$target
#id2info->{$id}->{'start'}=$start
#id2info->{$id}->{'end'}=$end
###################################################
sub getProteinFromBlastx_old{
	my ($blastOutFile)=@_;

	my $id2seq;
	my $id2info;
	my $oneseq="";
	open(IN,$blastOutFile)or die"$!";
	while(<IN>){
		chomp;
		if(/^Query\=(.+)$/){
			my $newid=$1;
			if(length($oneseq) > 10){
				$id2seq->{$id}=$oneseq;
			}
			$id=$newid;
			$id=~s/^\s{1,}//;
			$id=~s/\s{1,}$//;
			$oneseq="";
		}elsif(/^>(.+)$/){
			my $target=$1;
			$target=~s/^\s{1,}//;
			$target=~s/\s{1,}$//;
			$id2info->{$id}->{'target'}=$target;
		}elsif(/^Query\:/){
			my @line=split(/\s{1,}/);
			if(exists $id2info->{$id}->{'start'}){  #如果存在一个query对应的target中有两个区域都匹配，则需要考虑这个问题
				my $startOld=$id2info->{$id}->{'start'}; 
				if($startOld > $line[1]){
					next;
				}
			}else{
				$id2info->{$id}->{'start'}=$line[1];
			}
			$id2info->{$id}->{'end'}=$line[3];
			$oneseq.=$line[2];
			
		}
	}
	close IN or die"$!";

	if(length($oneseq) > 10){#处理最后一条记录
		$id2seq->{$id}=$oneseq;
	}

	return ($id2seq,$id2info);
}


################################################################################
#This is to read the fasta sequences;
#The input is the file name;
#The output is the sequence variable,with id linked to one sequence;
################################################################################
sub readSeq{
    my($infile)=@_;
    open(IN,$infile)or die"$!";
    my $seqRef;
    my $id;
    while(<IN>){
       chomp;
       next if(/^$/);
       if(/^>(.+)$/){
           $id=$1;
           next;
        }
       my @line=split(//);
       push(@{$seqRef->{$id}},@line);
     }
    close IN or die"$!";

    return($seqRef);
}


###############################################################################
#从blast的输出结果中，给出每个查询序列能够在数据库中找到的最相似的target的score
#输入为blast的输出文件路径，输出为每个序列的ID对应的score，有些ID可能没有找到hit
###############################################################################

sub getScoreFromBlast{
	my ($blastOutfile)=@_;

	open(IN,$blastOutfile)or die"$!";
	my $id2score;
	my $id;
	while(<IN>){
		chomp;
		next if(/^$/);
		if(/Query\=(.+?)$/){
			$id=$1;
			$id=~s/^\s{1,}//;
			$id=~s/\s{1,}$//;
		}elsif(/Sequences producing significant alignments/){
			while(<IN>){
				chomp;
				next if(/^$/);	
				last if(/^>/);
				if(/^[a-z0-9]/i){
					my @line=split(/\s{1,}/);
					my $score=$line[1];
					$id2score->{$id}=$score;
				}
			}
		}
	}
	close IN or die"$!";

	return $id2score;
}


########################################################################################
#此函数是给定序列的ID和序列，输出该序列在6种编码框中可能编码的所有长度大于50的蛋白序列
#输入为序列的ID和序列
#输出为肽段序列
#########################################################################################

sub translate6{
	my ($minLen,@seq)=@_;
	my $temp=join("",@seq);
	$temp=~s/^\-{1,}//;
	$temp=~s/\-{1,}$//;
	@seq=split(//,$temp);

	my @ORF;
	for(my $step=1;$step<4;$step++){
		my $start=$step;
		my $orf="";
		my $tag=0;
		for(my $i=$start-1;$i<=scalar @seq-3;$i+=3){
			my $codon=$seq[$i].$seq[$i+1].$seq[$i+2];
			$codon=uc($codon);
			next if($codon=~/---/);
			if($codon=~/TAA/ or $codon=~/TAG/ or $codon=~/TGA/){
				if($orf ne ''){
					push(@ORF,$orf);
				}
				$orf="";
			}else{
				$orf.=$codon;
				if($i >= (scalar @seq -5)){
					if($orf ne ''){
						push(@ORF,$orf);
					}
					$orf="";
				}
			}
		}
	}
	
	@seq=reverse @seq;
	$temp=join("",@seq);
	$temp=~tr/[ATCG]/[TAGC]/;
	@seq=split(//,$temp);  #这里是得到互补链
	for(my $step=1;$step<4;$step++){
		my $start=$step;
		my $orf="";
		my $tag=0;
		for(my $i=$start-1;$i<=scalar @seq-3;$i+=3){
			my $codon=$seq[$i].$seq[$i+1].$seq[$i+2];
			$codon=uc($codon);
			next if($codon=~/---/);
			if($codon=~/TAA/ or $codon=~/TAG/ or $codon=~/TGA/){
				push(@ORF,$orf);
				$orf="";
			}else{
				$orf.=$codon;
				if($i >= (scalar @seq -5)){
					push(@ORF,$orf);
					$orf="";
				}
			}
		}
	}

	my @peptide;
	foreach my$orf(@ORF){
		if(length($orf) < 3*$minLen){
			next;
		}
		my $onePep=DNA2protein($orf,0);
		$onePep=~s/^\-{1,}//;
		$onePep=~s/\-{1,}$//;
		push(@peptide,$onePep);
	}

	return @peptide;
}



################################################################################################
# A subroutine to translate a DNA 3-character codon to an amino acid

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'STOP',    # Stop
    'TAG' => 'STOP',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'STOP',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
     );

    if(exists $genetic_code{$codon}) {
          return $genetic_code{$codon};
    }elsif($codon =~/-/){
          return '-';
    }elsif(length($codon)!=3){
          print "$codon length is wrong\n";
          return 'error';
    }else{
          print  "Bad codon \"$codon\"!!\n";
          return 'X';  
    }
}


#################################################################################################
#This is to translate a DNA to protein.

sub DNA2protein{
     my ($dna,$starting_position)=@_;
     if($starting_position eq ''){
         $starting_position=0;
      }
     my $protein = '';

     for(my $i=$starting_position; $i < (length($dna)- 2) ; $i += 3) {                                             
         my $aa=codon2aa(substr($dna,$i,3));
         if($aa eq 'STOP'){
                last;
         }elsif($aa eq 'error'){
                last;
         }else {
               $protein .=$aa;
          }
    }
    return $protein;
}



###################################################End################################
1;
