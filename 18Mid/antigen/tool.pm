









##############################################################################
#该脚本是给定病毒之间的aa差异和其他一些变量值，给出这些病毒之间的ratio值。
#其输入包括下面几个：
#1 $dataset   => $dataset->[$#indexes]=\@compare
#2 $indexes  => @{$indexes}   #viral pairs
#3 $groups   => @{$groups}  #antigenic epitopes
#4 $site_groups => $site_groups->{$site}=$group   #the epitope each site belongs to
#5 $chars  => @{$chars}  #the aa index
#6 $character => $character->{$chars}->{$aa}=value  #specify the value of aa in the aa index database
#7 ...
#输出为每对病毒之间的ratio值，格式为: A\tB\tratio
#################################################################################

sub AV{
	my ($dataset,$indexes,$groups,$site_groups,$chars,$character,$receptor,$threshold,$prior,$distribution,$glyc)=@_;
	my @indexes=@{$indexes};

	#	Matrix
	print "\t...Matrix...\n";
	my $data;
	my $index=0;

	my $groups_data; #$groups->[group index]->[case index]=changed sites number
	$groups_data=Groups($dataset,$groups,$site_groups);
	for(my $i=0;$i<=$#{$groups_data};$i++){
		@{$data->[$index]}=@{$groups_data->[$i]};
		$index++;
	}

	my $chars_data; #$chars->[char index]->[case index]=changed aa character value
	$chars_data=Chars($dataset,$chars,$character);
	for(my $i=0;$i<=$#{$chars_data};$i++){
		@{$data->[$index]}=@{$chars_data->[$i]};
		$index++;
	}

	my @receptor;
	@receptor=Receptor($dataset,$receptor);
	@{$data->[$index]}=@receptor;
	$index++;

	my @glyc;
	@glyc=Glyc(\@indexes,$glyc);
	@{$data->[$index]}=@glyc;

	#	Prediction
	my @predRatio;  # 储存单次预测的结果
	for(my $i=0;$i<=$#indexes;$i++){
		my @input;
		for(my $j=0;$j<=$#{$data};$j++){
			if($data->[$j]->[$i] > $threshold->[$j]){
				push(@input,1);
			}else{
				push(@input,0);
			}
		}
		my $ratio = Prediction($distribution,@input);
		$ratio*=$prior;
		push(@predRatio,$indexes[$i]."\t".$ratio);
	}

	return (@predRatio);
}



########################################################################
#此程序是预测序列中的糖基化位点(Asn-Xaa-Ser/Thr,Xaa为除了Pro之外的aa）。
#输入为序列，为字符串的形式；
#输出为糖基化位点的起始位置,为数组；
#######################################################################
sub ngly{
	my ($seq)=@_;
	my @pos;
	while($seq=~/N[^P][ST]/g){
		my $pos=pos $seq;
		pos($seq)=$pos-2;  #为了对所有的子字符进行匹配
		my $new=$pos-2;
		if(scalar @pos > 0){
			if($pos[-1] == $new-1){  #如果之前的位置也是糖基化，那么之前的位置用现在的代替
				$pos[-1]=$new;
			}else{
				push(@pos,$new);
			}
		}else{
			push(@pos,$new);
		}
	}
	return @pos;
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
		s/^[\s\r\n\t]{1,}//;
		s/[\r\t\s\n]{1,}$//;
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



##################################################################################################
#此函数是把用户的序列整理成与参考序列相似的形式
#输入为需要处理的序列文件，输出文件名，参考序列，临时文件夹
#注意：参考序列的名字最好跟待处理序列的名字有差异！！！
#没有输出
#这里需要注意的是，如果参考序列本身在某些位点有gap，可能会导致最终处理的序列在这些位点上缺失。因此，这里采用了
#一个简单的方法，即那些通过比对引进的位点在参考序列的位点中应该绝大部分都是gap，这里要求gap比例超过90%。
##################################################################################################
sub align2referSeq{
	my ($seqForAlignFile,$seqAlignFile,$referSeqFile,$tempDir)=@_;
	
	my $combinedFile=$tempDir."combineSeq";
	my $combinedFile_align=$combinedFile."_align";
	`cat $seqForAlignFile $referSeqFile > $combinedFile`;
	`mafft --retree 1 $combinedFile > $combinedFile_align`;

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
		my $largestRatio=$siteRef_standard->{$site}->{$key[0]}/(scalar @referID);
        if($key[0] eq '-' && $largestRatio > 0.95){  #如果标准序列的某个位点是空格，那么该位置应该去掉，除非一些特殊的位点
			next;
		}else{
			push(@remainPos,$site);
		}
    }

	#finally, get the aligned sample sequences
	open(OUT,">$seqAlignFile")or die"$!";
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







#################################################################################################
1;

