############################################################################################
#	Purpose : 0) Cross Validation 1) data manage for each attributes 2)discrete every attribute, 3) construct Naive Bayesian using Bernoulli models and uniform proior distribution, 4) prediction
#	<For 2008.6.2>
##

##
#	sub function : Model
#	Purpose : Read Model From Files
#	Input : $infile
#	Output : @groups,%groups,@chars,$character,%receptor,@threshold,$distribution
################################################################################################
sub Model
{
	use strict;
	my ($file) = @_;
	my @groups;
	my %groups;
	my @chars;
	my $character;
	my %receptor;
	my @threshold;
	my $distribution;
	my $prior;
	
	open(IN,"<$file")||die "$!";
	while(<IN>){
		chomp;
		next if(/^$/);
		if(/^\<Attribute:Group\>$/){
#			print "\tAttribute:Group...\n";
			while(<IN>){
				last if(/^\/\/$/);
				s/^\s+//;
				my @lineArray=split(/\s+/);
				my $group=shift @lineArray;
				push(@groups,$group);
				foreach my $one(@lineArray){
					$groups{$one}=$group;
				}
			}
		}elsif(/^\<Attribute:Char\>$/){
#			print "\tAttribute:Char...\n";
			$_=<IN>;
			s/^\s+//;
			my @aa=split(/\s+/);
			shift @aa;
			
			while(<IN>){
				last if(/^\/\/$/);
				s/^\s+//;
				my @lineArray=split(/\s+/);
				my $char=shift @lineArray;
				die "Bad line formate in Char context : $_\n" if($#lineArray != $#aa);
				push(@chars,$char);
				for(my $i=0;$i<=$#lineArray;$i++){
					$character->{$char}->{$aa[$i]}=$lineArray[$i];
				}
			}
		}elsif(/^\<Attribute:Receptor\>$/){
#			print "\tAttribute:Receptor...\n";
			while(<IN>){
				last if(/^\/\/$/);
				s/^\s+//;
				my @lineArray=split(/\s+/);
				die "Bad line formate in Receptor context : $_\n" if($#lineArray != 1);
				$receptor{$lineArray[0]}=$lineArray[1];
			}
		}elsif(/^\<Parameter\>$/){
#			print "\tParameter...\n";
			while(<IN>){
				last if(/^\/\/$/);
				s/^\s+//;
				my @lineArray=split(/\s+/);
				die "Bad line formate in Parameter context : $_\n" if($#lineArray != 1);
				$threshold[$lineArray[0]]=$lineArray[1];
			}
		}elsif(/^\<Model\>$/){
#			print "\tModel...\n";
			$_=<IN>;
			s/^\s+//;
			my @lineArray=split(/\s+/);
			die "Bad line formate for prior : $_\n" if($#lineArray != 1);
			$prior=$lineArray[1];
			while(<IN>){
				last if(/^\/\/$/);
				s/^\s+//;
				my @lineArray=split(/\s+/);
				die "Bad line formate in Model context : $_\n" if($#lineArray != 4);
				$distribution->[$lineArray[0]]->[0]->[0]=$lineArray[1];
				$distribution->[$lineArray[0]]->[0]->[1]=$lineArray[2];
				$distribution->[$lineArray[0]]->[1]->[0]=$lineArray[3];
				$distribution->[$lineArray[0]]->[1]->[1]=$lineArray[4];
			}
		}
	}
	close(IN)||die "$!";
	return(\@groups,\%groups,\@chars,$character,\%receptor,\@threshold,$prior,$distribution);
}

####################################################################################
#	sub function : CV
#	Purpose : Cross Validation
#	Input : $dataset,$cases(index number),$fold(number)
#	Output : $accuracy
###################################################################################
sub CV
{
	use strict;
	my ($dataset,$fold,@marks) = @_;
	my $accuracy=0;
	my @indexes=(0..$#marks);
	
	die "Fold number bigger than number of cases ( ".scalar(@marks)." )!\n" if($fold > @marks);
	
	my @rand;
	while($#indexes > 0){
		my $rand=int(rand(scalar @indexes));
		push(@rand,$indexes[$rand]);
		splice(@indexes,$rand,1);
	}
	push(@rand,@indexes);

	my $interval=scalar(@rand)/$fold;
	
	my $cv=[];
	for(my $i=0;$i<$fold;$i++){
		if($i+1 == $fold){
			@{ $cv->[$i] } = @rand[$i*$interval..$#rand];
		}else{
			@{ $cv->[$i] } = @rand[$i*$interval..(($i+1)*$interval-1)];
		}
	}

	print "\t$fold fold cross validation ...\n";
	for(my $i=0;$i<=$#{$cv};$i++){
		print "\t\t".($i+1)." ...\t";
		my $train=[];
		my @train_marks;
		my $test=[];
		my @test_marks=@marks[@{$cv->[$i]}];

		print " [1] Dataset";
		for(my $k=0;$k<=$#{$cv};$k++){
			if($k==$i){
				for(my $j=0;$j<=$#{$dataset};$j++){
					@{$test->[$j]}=@{ $dataset->[$j] }[ @{$cv->[$i]} ];
				}
			}else{
				for(my $n=0;$n<=$#{$cv->[$k]};$n++){
					for(my $j=0;$j<=$#{$dataset};$j++){
						push(@{$train->[$j]},$dataset->[$j]->[ $cv->[$k]->[$n] ]);
					}
					push(@train_marks,@marks[ $cv->[$k]->[$n] ]);
				}
			}
		}

		##	Discretion
		print " [2] Discretion";
		my @discretion = Discretion($train,@train_marks);

		##	Learning
		print " [3] Learning";
		my $distribution = NaiveBayesian($train,\@train_marks,\@discretion);

		##	Testing
		print " [4] Testing\n";
		for(my $k=0;$k<=$#test_marks;$k++){
			my @input;
			for(my $j=0;$j<=$#{$test};$j++){
				if($test->[$j]->[$k] > $discretion[$j]){
					push(@input,1);
				}else{
					push(@input,0);
				}
			}
			my $ratio = Prediction($distribution,@input);

			my $change=0;
			foreach(@train_marks){
				$change+=$_;
			}
			$ratio*=(1+$change)/(1+scalar(@train_marks)-$change);

			my $mark=0;
			$mark=1 if($ratio > 1);
			$accuracy++ if($mark == $test_marks[$k]);
		}
	}
	$accuracy/=scalar(@marks);
	return $accuracy;
}


####################################################################################
#与CV（)的差别在于会给出准确率，敏感性和特异性
#	sub function : CV2
#	Purpose : Cross Validation
#	Input : $dataset,$cases(index number),$fold(number)
#	Output : $accuracy,$sensitivity,$specificity and the prediction values;
###################################################################################
sub CV2
{
	use strict;
	my ($dataset,$fold,@marks) = @_;
	my @indexes=(0..$#marks);
	
	die "Fold number bigger than number of cases ( ".scalar(@marks)." )!\n" if($fold > @marks);
	
	my @rand;
	while($#indexes > 0){
		my $rand=int(rand(scalar @indexes));
		push(@rand,$indexes[$rand]);
		splice(@indexes,$rand,1);
	}
	push(@rand,@indexes);

	my $interval=scalar(@rand)/$fold;
	
	my $cv=[];
	for(my $i=0;$i<$fold;$i++){
		if($i+1 == $fold){
			@{ $cv->[$i] } = @rand[$i*$interval..$#rand];
		}else{
			@{ $cv->[$i] } = @rand[$i*$interval..(($i+1)*$interval-1)];
		}
	}

	print "\t$fold fold cross validation ...\n";
	my (@prediction,@predict,@obs);  #@prediction store the values and @predict store the '0' or '1';

	for(my $i=0;$i<=$#{$cv};$i++){
		print "\t\t".($i+1)." ...\t";
		my $train=[];
		my @train_marks;
		my $test=[];
		my @test_marks=@marks[@{$cv->[$i]}];

		print " [1] Dataset";
		for(my $k=0;$k<=$#{$cv};$k++){
			if($k==$i){
				for(my $j=0;$j<=$#{$dataset};$j++){
					@{$test->[$j]}=@{ $dataset->[$j] }[ @{$cv->[$i]} ];
				}
			}else{
				for(my $n=0;$n<=$#{$cv->[$k]};$n++){
					for(my $j=0;$j<=$#{$dataset};$j++){
						push(@{$train->[$j]},$dataset->[$j]->[ $cv->[$k]->[$n] ]);
					}
					push(@train_marks,@marks[ $cv->[$k]->[$n] ]);
				}
			}
		}

		##	Discretion
		print " [2] Discretion";
		my @discretion = Discretion($train,@train_marks);

		##	Learning
		print " [3] Learning";
		my $distribution = NaiveBayesian($train,\@train_marks,\@discretion);

		##	Testing
		print " [4] Testing\n";
		for(my $k=0;$k<=$#test_marks;$k++){
			my @input;
			for(my $j=0;$j<=$#{$test};$j++){
				if($test->[$j]->[$k] > $discretion[$j]){
					push(@input,1);
				}else{
					push(@input,0);
				}
			}
			my $ratio = Prediction_test($distribution,@input);
#			my $ratio = Prediction($distribution,@input);

			my $change=0;
			foreach(@train_marks){
				$change+=$_;
			}
			$ratio*=(1+$change)/(1+scalar(@train_marks)-$change);
			push(@prediction,$ratio);

			my $mark=0;
			$mark=1 if($ratio > 1);
			push(@predict,$mark);
			push(@obs,$test_marks[$k]);
		}
	}
	my ($accuracy,$sensitivity,$specificity)=evaluate2(\@obs,\@predict);
	return ($accuracy,$sensitivity,$specificity,\@prediction,\@obs);
}

##############################################################################################
#forROC()   #给定实际值和预测值，得到随着cutoff的变化而变化的true positive rate和false positive rate;
# Input: Two arrays, one is the observed and the other is the predicted data;
# Output: Two arrays, one is the TPR and the other is the FPR;
##############################################################################################
sub forROC{
	my ($obs,$pre)=@_;
	my @obs=@{$obs};
	my @pre=@{$pre};
	my $cutoffRef;
	foreach my$one(@{$pre}){
		$cutoffRef->{$one}++;
	}

	my (@TPR,@FPR);
	my @cutoff=sort {$a<=>$b}keys %{$cutoffRef};
	for(my $i=0;$i<scalar @cutoff;$i++){
		my @pre_discre;
		for(my $j=0;$j<scalar @pre;$j++){
			if($pre[$j] >= $cutoff[$i]){
				push(@pre_discre,1);
			}else{
				push(@pre_discre,0);
			}
		}
		my ($agree,$sen,$spe)=evaluate2(\@obs,\@pre_discre);
		push(@TPR,$sen);
		push(@FPR,1-$spe);
	}
	return(\@TPR,\@FPR);
}


#################################################################################################
#evaluate2()     #calculate the sensitivity and specificity of the prediction result;
#The input is the observed and predicted data.Attention:The observed should be in first.
#The output is the agreement rate,sensitivity and specificity.
#该程序和evaluate的不同在于，输入的是二元数据，即0和1.如果预测的和实际的都是1或者都是0，则预测正确。
#################################################################################################
sub evaluate2{
    my($obs,$pre)=@_;
    my @obs=@$obs;
    my @pre=@$pre;
    if(scalar @obs != scalar @pre){
       die"Error!The number of observations is not equal to that of predicted\n";
     }
    if(scalar @obs ==0){
       die"No observations!The array is empty!\n";
     }
    my($posTotal,$negTotal,$agree,$sen,$spe)=(0,0,0,0,0);
    for(my $i=0;$i<scalar @obs;$i++){
       die unless($obs[$i]=~/[01]/);
       die unless($pre[$i]=~/[01]/);
       if($obs[$i] == $pre[$i]){
          $agree++;
         }
       if($obs[$i] == 1){
          $posTotal++;
          if($pre[$i] ==1){
             $sen++;
            }
         }
       if($obs[$i] ==0){
          $negTotal++;
          if($pre[$i] ==0){
             $spe++;
            }
         }
     }
    return($agree/(scalar @obs),$sen/$posTotal,$spe/$negTotal);
}
 

#############################################################################################
#	sub Function : Groups
#	Purpose : Storing the Group Feature Information for Later Bayesian Network Analysis
#	Input : $dataset($dataset->[case index]=\@changes),@groups(List of groups),%groups($groups{site}=group)
#	Output : $data($data->[group index]->[case index]=site number changed)
##############################################################################################
sub Groups
{
	use strict;
	my ($dataset,$attribute,$assignment) = @_;
	my $data;

	for(my $i=0;$i<=$#{$dataset};$i++){
		my %distribution;
		my @mutations=@{ $dataset->[$i] };
		@mutations=map /^\D+(\d+)\D+$/,@mutations;
		foreach my $one(@mutations){
			$distribution{$assignment->{$one}}++ if(exists $assignment->{$one});
		}
		for(my $j=0;$j<=$#{$attribute};$j++){
			if(exists $distribution{$attribute->[$j]}){
				$data->[$j]->[$i] = $distribution{ $attribute->[$j] };
			}else{
				$data->[$j]->[$i] = 0;
			}
		}
	}

	return $data;
}

##################################################################################################
#	sub Function : Chars
#	Purpose : Storing the aa character Feature Information for Later Bayesian Network Analysis
#	Input : $dataset($dataset->[case index]=\@changes),@chars(List of chars),%character($character->{index}->{aa}=value)
#	Output : $data($data->[char index]->[case index]=aa character changed value)
###################################################################################################
sub Chars
{
	use strict;
	my ($dataset,$attribute,$assignment) = @_;
	my $data;

	for(my $i=0;$i<=$#{$dataset};$i++){
		my $changes; #$changes->{attribute}=\@changed values
		foreach my $one(@{ $dataset->[$i] }){
			my ($a1,$site,$a2) = ($one =~ /^(\w)(\d+)(\w)$/);
			foreach my $two(@{$attribute}){
				push(@{$changes->{$two}},abs( $assignment->{$two}->{$a1} - $assignment->{$two}->{$a2} ) );
			}
		}
		for(my $j=0;$j<=$#{$attribute};$j++){
			if($changes->{$attribute->[$j]} eq ""){   #对于每个特征，选择改变最大的三个位点的改变量的均值作为该特征的改变；
				push(@{$data->[$j]},0);
				next;
			}
			my @changes=reverse sort {$a <=> $b} @{$changes->{$attribute->[$j]}};
			my $distance=0;
			my $maxIndex=2;
			if($#changes == -1){
				push(@{ $data->[$j] },0);
				next;
			}

			$maxIndex=$#changes if($#changes < $maxIndex);
			foreach(@changes[0..$maxIndex]){
				$distance+=$_;
			}
			$distance/=($maxIndex+1);
			push(@{ $data->[$j] },$distance);
		}
	}

	return $data;
}

#############################################################################################
#	sub Function : Receptor
#	Purpose : Storing the receptor binding influence Feature Information for Later Bayesian Network Analysis
#	Input : $dataset($dataset->[case index]=\@changes),%features($features{$site}=receptor binding influence value)
#	Output : @data(for each cases)
############################################################################################
sub Receptor
{
	use strict;
	my ($dataset,$assign) = @_;
	my @data;
	my %assignment=%{$assign};

	for(my $i=0;$i<=$#{$dataset};$i++){
		my @changes;
		my @mutations=@{ $dataset->[$i] };
		if($#mutations == -1){
			push(@data,0);
			next;
		}
		@mutations=map /^\D+(\d+)\D+$/,@mutations;
		foreach my $one(@mutations){
			push(@changes,$assignment{$one});
		}
		@changes = reverse sort {$a <=> $b} @changes;
		my $distance=0;
		my $maxIndex=2;
		if($#changes == -1){
			push(@data,0);
			next;
		}
		$maxIndex=$#changes if($#changes < $maxIndex);
		foreach(@changes[0..$maxIndex]){
			$distance+=$_;
		}
		$distance/=($maxIndex+1);
		push(@data,$distance);
	}

	return @data;

}

###########################################################################################
#	sub Function : Glyc
#	Purpose : Storing the glycosylation site Feature Information for Later Bayesian Network Analysis
#	Input : $glyc($glyc->{sequnce}=\@glycosylation sites),@cases(List of index pair)
#	Output : @data(for each cases)
###########################################################################################
sub Glyc
{
	use strict;
	my ($dataset,$assignment) = @_;
	my @data;

	foreach(@{$dataset}){
		my @indexPair=split(/\t/);
		die "Bad index pair formate (index1".'\t'."index2)!\n" if($#indexPair != 1);
		my %seen;
		foreach(@{ $assignment->{ $indexPair[0] } }){
			$seen{$_}=1;
		}
		my @same=grep $seen{$_},@{ $assignment->{ $indexPair[1] } };
		push(@data,scalar(@{$assignment->{$indexPair[0] }}) + scalar(@{$assignment->{ $indexPair[1] }}) - 2*scalar(@same) );
	}

	return @data;
}

############################################################################################
#	sub Function : Dataset
#	Purpose : Store the training and testing dataset in memory
#	Input : @train/@test,@groups,%groups
#	Output : @antigenicity,$data->[group index]->[case index]=$value
###########################################################################################
sub Dataset
{
	use strict;
	my ($dataset,$attribute,$assignment) = @_;
	my @antigenicity;
	my $data=[];

	for(my $i=0;$i<=$#{$dataset};$i++){
		my %distribution;
		my @mutations=split(/\t/,$dataset->[$i]);
		push(@antigenicity,shift @mutations);
		@mutations=map /^\D+(\d+)\D+$/,@mutations;
		foreach my $one(@mutations){
			$distribution{$assignment->{$one}}++ if($assignment->{$one} ne "");
		}
		for(my $j=0;$j<=$#{$attribute};$j++){
			if($distribution{$attribute->[$j]} ne ""){
				$data->[$j]->[$i] = $distribution{ $attribute->[$j] };
			}else{
				$data->[$j]->[$i] = 0;
			}
		}
	}
	return ($data,@antigenicity);
}

############################################################################################
#	sub Function : Discretion
#	Purpose : Discrete every attribute
#	Input : $data,@antigenicity
#	Output : @discretion($discretion[attribute index]=threshold)
############################################################################################
sub Discretion
{
	use strict;
	my ($dataset,@label) = @_;
	my @discretion;

	for(my $i=0;$i<=$#{$dataset};$i++){

		my $threshold;
		my $value=0;

		my @changes=@{$dataset->[$i]};
		my @sortArray;
		foreach my $one(@changes){
			push(@sortArray,$one) if(!grep $_ == $one,@sortArray);
		}
		@sortArray=sort {$a <=> $b} @sortArray;
		if($#sortArray == 0){
			$discretion[$i]=$sortArray[0];
			next;
		}

		for(my $j=0;$j<=$#sortArray;$j++){
			my $N=[];
			$N->[0]->[0]=$N->[0]->[1]=$N->[1]->[0]=$N->[1]->[1]=0;
			for(my $k=0;$k<=$#label;$k++){
				if( ($label[$k] == 0) && ($changes[$k] <= $sortArray[$j]) ){
					$N->[0]->[0]++;
				}elsif( ($label[$k] == 0) && ($changes[$k] > $sortArray[$j]) ){
					$N->[0]->[1]++;
				}elsif( ($label[$k] == 1) && ($changes[$k] <= $sortArray[$j]) ){
					$N->[1]->[0]++;
				}else{
					$N->[1]->[1]++;
				}
			}
			my $valueTem = ChiSquare($N);
			if($valueTem > $value){
				$threshold=$sortArray[$j];
				$value=$valueTem;
			}
		}
		if($threshold eq ""){
			$discretion[$i]=$sortArray[0];
		}else{
			$discretion[$i]=$threshold;
		}
	}
	return @discretion;
}

###################################################################################################
#	sub Function : NaiveBayesian
#	Purpose : Parameter Learning
#	Input : $data,@antigenicity,@discretion
#	Output : $distribution->[attribute index]->[0/1 for attribute]->[0/1 for antigenicity]=$count
#注意：这里输出的是属性值在前，抗原变化在后的数据分布
#################################################################################################
sub NaiveBayesian
{
	use strict;
	my ($dataset,$label,$discretion) = @_;
	my $distribution=[];
	
	for(my $i=0;$i<=$#{$dataset};$i++){
		$distribution->[$i]->[0]->[0]=0;
		$distribution->[$i]->[0]->[1]=0;
		$distribution->[$i]->[1]->[0]=0;
		$distribution->[$i]->[1]->[1]=0;
		for(my $j=0;$j<=$#{$dataset->[$i]};$j++){
			my $index;
			if($dataset->[$i]->[$j] > $discretion->[$i]){
				$index=1;
			}else{
				$index=0;
			}
			$distribution->[$i]->[$index]->[ $label->[$j] ]++;
		}
	}
	return $distribution;
}

##############################################################################
#	sub Function : Prediction
#	Purpose : Make Naive Bayesian Using Learned Parameter
#	Input : @new,$distribution
#	Output : $probability
###############################################################################
sub Prediction
{
	use strict;
	my ($distribution,@new) = @_;
	my $result=1;
	
	for(my $i=0;$i<=$#{$distribution};$i++){
		my ($n0,$n1);
		$n0=$distribution->[$i]->[0]->[0] + $distribution->[$i]->[1]->[0];
		$n1=$distribution->[$i]->[0]->[1] + $distribution->[$i]->[1]->[1];
		my $resultTem=( (1+ $distribution->[$i]->[ $new[$i] ]->[1])/(2 + $n1) )/( (1+ $distribution->[$i]->[ $new[$i] ]->[0])/(2 + $n0) );
		$result*=$resultTem;
	}
	return $result;
}

##############################################################################
#此程序是为了测试各个属性的效果
###############################################################################
sub Prediction_test
{
	use strict;
	my ($distribution,@new) = @_;
	my $result=1;
	
	for(my $i=0;$i<=$#{$distribution};$i++){
		next if($i ==0 or $i==2 or $i==3 or $i==10);
		my ($n0,$n1);
		$n0=$distribution->[$i]->[0]->[0] + $distribution->[$i]->[1]->[0];
		$n1=$distribution->[$i]->[0]->[1] + $distribution->[$i]->[1]->[1];
		my $resultTem=( (1+ $distribution->[$i]->[ $new[$i] ]->[1])/(2 + $n1) )/( (1+ $distribution->[$i]->[ $new[$i] ]->[0])/(2 + $n0) );
		$result*=$resultTem;
	}
	return $result;
}


##############################################################################
#此程序是为了测试当只用某些属性的时候模型的效果
###############################################################################
sub Prediction_test2
{
	use strict;
	my ($distribution,$removeRef,@new) = @_;
	my $result=1;
	
	for(my $i=0;$i<=$#{$distribution};$i++){
		next if(exists $removeRef->{$i});
		my ($n0,$n1);
		$n0=$distribution->[$i]->[0]->[0] + $distribution->[$i]->[1]->[0];
		$n1=$distribution->[$i]->[0]->[1] + $distribution->[$i]->[1]->[1];
		my $resultTem=( (1+ $distribution->[$i]->[ $new[$i] ]->[1])/(2 + $n1) )/( (1+ $distribution->[$i]->[ $new[$i] ]->[0])/(2 + $n0) );
		$result*=$resultTem;
	}
	return $result;
}

###############################################################################
#	sub function : ChiSquare
#	Purpose : Genearte chi-sequare test value
#	Input : 2x2 contigency table
#	Output : value
###############################################################################
sub ChiSquare
{
	use strict;
	my ($table) = @_;
	my $value=0;

	my ($N,@Ni,@Nj);
	$N=$table->[0]->[0]+$table->[0]->[1]+$table->[1]->[0]+$table->[1]->[1];
	@Ni=($table->[0]->[0]+$table->[0]->[1],$table->[1]->[0]+$table->[1]->[1]);
	@Nj=($table->[0]->[0]+$table->[1]->[0],$table->[1]->[1]+$table->[0]->[1]);

	for(my $i=0;$i<=$#{$table};$i++){
		for(my $j=0;$j<=$#{$table->[$i]};$j++){
			if($Ni[$i]*$Nj[$j] == 0){
				$value+=0;
			}else{
				$value+=($table->[$i]->[$j] - $Ni[$i]*$Nj[$j]/$N)*($table->[$i]->[$j] - $Ni[$i]*$Nj[$j]/$N)/($Ni[$i]*$Nj[$j]/$N);
			}
		}
	}
	return $value;
}
1;
