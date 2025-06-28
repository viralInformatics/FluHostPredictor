#!/usr/bin/perl

####################################################################################
# Author:  KM Amada (kmamada@ifrec.osaka-u.ac.jp)
#
# Ver. Date      Changelog
####################################################################################
# 1.0  11.01.13  Initial release
#
# **Skipped version 2 to standardise version numbers to seekquencer.pl script**
#
# 3.0  04.24.14  Added split option -mod <mafftash-split> for output
#                Uses seekquencer_v3 backend
#
# 4.0  05.12.14  Added new options: -run <thread|normal> -trd <count> -noin
#                Sets -seqa fast in seekquencer.pl
#                Uses seekquencer_v4 backend
#
# 4.1  05.19.14  Added a check on running REST requests before proceeding
#                to avoid server load problems
#
####################################################################################

use strict;
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use LWP::Simple;
use LWP::UserAgent;

# to prevent error: Header line too long (limit is 8192)
use LWP::Protocol::http;
push(@LWP::Protocol::http::EXTRA_SOCK_OPTS, MaxLineLength => 0);



my $BASEURL = "http://sysimm.ifrec.osaka-u.ac.jp/seekquencer/REST/service.cgi/premafft";
my ( $INPUTFILE, $IDLISTFILE, $SEQFASTAFILE, $OUTPUTFILE, $SEQFLAG, $STRFLAG, $EVALFLAG, $NOINFLAG );
my $OUTTYPE = "mafftash";
my $SEQLIMIT = 100;

my $RUNMODE = "normal";    # thread|normal
my $THREADCOUNT = 3;


GetOptions
(
    'inp=s' => \$INPUTFILE,
    'idf=s' => \$IDLISTFILE,
    'seqf=s'=> \$SEQFASTAFILE,
    'out=s' => \$OUTPUTFILE,
    'str'   => \$STRFLAG,
    'seq'   => \$SEQFLAG,
    'lim=i' => \$SEQLIMIT,
    'pre'   => \$EVALFLAG,
    'noin'  => \$NOINFLAG,
    'mod=s' => \$OUTTYPE,
    'run=s' => \$RUNMODE,
    'trd=i' => \$THREADCOUNT,


);

print STDERR "[Seekquencer-premafft 4.1]\n";

# set temp directory
my $TMP = "/tmp/seekpremafft$$";
make_path($TMP) unless -d $TMP;



######
# validation
help("Required parameter: define input as '-inp' or '-idf' or '-seqf'") if ( !defined $INPUTFILE && !defined $IDLISTFILE && !defined $SEQFASTAFILE );
help("'-inp' is already defined") if ( defined $INPUTFILE && (defined $IDLISTFILE || defined $SEQFASTAFILE) );
help("Input file $INPUTFILE does not exist (or filesize is 0)") if ( defined $INPUTFILE && (! -e $INPUTFILE || !-s $INPUTFILE) );
help("Input file $IDLISTFILE does not exist (or filesize is 0)") if ( defined $IDLISTFILE && (! -e $IDLISTFILE || !-s $IDLISTFILE) );
help("Input file $SEQFASTAFILE does not exist (or filesize is 0)") if ( defined $SEQFASTAFILE && (! -e $SEQFASTAFILE || !-s $SEQFASTAFILE) );
help("Required parameter: output file '-out'") unless ( defined $OUTPUTFILE );
help("Set either '-str' or '-seq' or dont set any at all") if ( defined $STRFLAG && defined $SEQFLAG );

help("Invalid value for '-mod <fasta|mafftash|mafftash-split>'") if ( $OUTTYPE ne "fasta" && $OUTTYPE ne "mafftash" && $OUTTYPE ne "mafftash-split" );
help("Invalid value for '-run <thread|normal>'") if ( $RUNMODE ne "thread" && $RUNMODE ne "normal" );
help("Invalid value for '-trd <count>'; count should be > 1") if ( $RUNMODE eq "thread" && $THREADCOUNT <= 1 && $THREADCOUNT > 10 );



######
# check existing requests
print STDERR "Checking server status...\n";

my $browser = LWP::UserAgent->new;
$browser->timeout(0);

# get: get results for mafftash job
my $jobsResponse = $browser->get("$BASEURL/jobs");

if ( $jobsResponse->is_success )
{
    my $runningJobs = parseJobQueryResponse($jobsResponse->content);
    print STDERR "Server Jobs: $runningJobs\n";

    bail("Max jobs reached. The server cannot process your request right now; try again later.", 0) if $runningJobs > 10;
}
else
{
    bail(sprintf("[%d] %s\n", $jobsResponse->code, parseError($jobsResponse->content)));
}


######
# make a temporary input if lists were provided
unless ( defined $INPUTFILE )
{
    $INPUTFILE = "$TMP/input.homemade";
    open INPF, ">$INPUTFILE" or bail("Error writing to input file.");

    if ( defined $IDLISTFILE )
    {
        open IDLIST, "<$IDLISTFILE" or bail("Error reading input file.");
        while( <IDLIST> )
    	{
        	chomp;
        	if ( /(\w{5})/ )
        	{
        	    print INPF ">PDBID\n$1\n";
        	}
        }
        close IDLIST;
    }


    if ( defined $SEQFASTAFILE )
    {
        open FASTA, "<$SEQFASTAFILE" or bail("Error reading input file.");
        while( <FASTA> )
    	{
        	chomp;
            print INPF "$_\n";
        }
        close FASTA;
    }

    close INPF;
}


######
# prepare parameters
print STDERR "Preparing parameters for service request...\n";

my @parameters = ();
push(@parameters, "fileinput" => ["$INPUTFILE"]);
push(@parameters, "out_type" => $OUTTYPE);

push(@parameters, "rest_flag" => "1");
push(@parameters, "cls_flag" => "1");
push(@parameters, "pre_flag" => "1") if defined $EVALFLAG;
push(@parameters, "noin_flag" => "1") if defined $NOINFLAG;

push(@parameters, "run_mode" => $RUNMODE);
push(@parameters, "thread_count" => $THREADCOUNT) if $RUNMODE eq "thread";


if ( defined $STRFLAG )
{
    push(@parameters, "str_flag" => "1");
    push(@parameters, "ash_flag" => "1");
}
elsif ( defined $SEQFLAG )
{
    push(@parameters, "seq_flag" => "1");
    push(@parameters, "seq_algorithm" => "fast");

    my $bfactor = $SEQLIMIT*10;
    push(@parameters, "seq_blastlimit" => $bfactor);
}
else
{
    push(@parameters, "str_flag" => "1");
    push(@parameters, "ash_flag" => "1");
    push(@parameters, "seq_flag" => "1");
    push(@parameters, "seq_algorithm" => "fast");

    my $bfactor = $SEQLIMIT*10;
    push(@parameters, "seq_blastlimit" => $bfactor);
}



######
# start rest service
print STDERR "Sending service request...\n";

# post: running a mafftash job
my $postResponse = $browser->post( $BASEURL, \@parameters, 'Content_Type' => 'form-data' );
bail(sprintf("[%d] %s\n", $postResponse->code, parseError($postResponse->content))) unless($postResponse->is_success);


# get response from post request
my ($status, $seekid) = parseResponse($postResponse->content);



my $MAXTRIES = 3;
my $STIMER = 4;
my $longtimer = 0;

print STDERR "Request sent! Waiting for response...[$seekid]\n";


# wait for results until it becomes available
while(1)
{
    $longtimer = $longtimer <= ($STIMER*3) ? $longtimer+$STIMER : $STIMER;
    sleep $longtimer;


    # get: get results for mafftash job
    my $getResponse = $browser->get("$BASEURL/$seekid");

    if ( $getResponse->is_success )
    {

        # get response from get request
        ($status, $seekid) = parseResponse($getResponse->content);
        next unless ( $status eq "done" );


        # if job is finished and ready
        print STDERR "Results found!\n";
        my $csfile = "$TMP/checksum.tar.gz";
        my $try1 = 1;


        while(1)
        {
            print STDERR "Fetching Results... [Trial $try1]\n";

            if ( is_success(getstore("$BASEURL/getmdlist/$seekid", $csfile)) && -e $csfile && -s $csfile )
            {
                # get response from get request
                my $checklist = extractchecksum($csfile);
                bail("Error retrieving list of compressed files!") unless ( scalar %$checklist > 0 );


                foreach my $id ( keys %$checklist )
                {
                    my $checkfile = "$TMP/$id";
                    my $checkid = $checklist->{$id};
                    my $try2 = 1;

                    while(1)
                    {
                        unlink $checkfile if -e $checkfile;

                        if ( is_success(getstore("$BASEURL/get/$seekid/$id", $checkfile)) && -e $checkfile && -s $checkfile )
                        {
                            my $hashid = getchecksum($checkfile);
                            #print STDERR "[hashid]$hashid [checkid]$checkid\n";

                            if ($hashid ne "" && $hashid ne $checkid )
                            {
                                unlink $checkfile if -e $checkfile;
                                bail("Error retrieving compressed file from server! [Checksum Failed]") if $try2 >= $MAXTRIES;
                                $try2++;
                                sleep $STIMER;
                            }
                            else
                            {
                                last;
                            }
                        }
                        else
                        {
                            bail("Error retrieving compressed file from server!") if $try2 >= $MAXTRIES;
                            $try2++;
                            sleep $STIMER;
                        }
                    }
                }

                last;
            }
            else
            {
                bail("Error retrieving list of compressed files from server!") if $try1 >= $MAXTRIES;
                $try1++;
                sleep $STIMER;
            }
        }

        last;

    }
    else
    {
        bail(sprintf("[%d] %s\n", $getResponse->code, parseError($getResponse->content)));
    }

}


# make sure outputs were generated
# decompress
print STDERR "Assembling final results...\n";
backticks("cat $TMP/archive.tar.gz* | tar -zxf - -C $TMP/");


if ( $OUTTYPE eq "mafftash-split" )
{
    my $xfile_str = "$TMP/$seekid.out.str";
    my $xfile_seq = "$TMP/$seekid.out.seq";
    bail("Error: Output file not found!") if ( !-e $xfile_str && !-e $xfile_seq );

    # just copy str output
    backticks("cp -f $xfile_str $OUTPUTFILE.str");

    # reprint seq output
    my $dataset = parseFasta($xfile_seq);
    my $seqctr = 0;

    open OUTFILE, ">$OUTPUTFILE.seq" or bail("Error writing to seq output file.");

    foreach my $ctr ( sort {$a <=> $b} keys %$dataset )
    {
        my $id  = $dataset->{$ctr}->{'id'};
        my $seq = $dataset->{$ctr}->{'sequence'};

        if ( $id !~ /^PDBID/ && $seqctr < $SEQLIMIT )
        {
            print OUTFILE ">$id\n$seq\n";
            $seqctr++;
        }
    }

    close OUTFILE;

}
elsif ( $OUTTYPE eq "mafftash" )
{
    my $xfile = "$TMP/$seekid.out";
    bail("Error: Output file not found!") unless -e $xfile;


    # reprint output
    my $dataset = parseFasta($xfile);
    my $seqctr = 0;

    open OUTFILE, ">$OUTPUTFILE" or bail("Error writing to output file.");

    foreach my $ctr ( sort {$a <=> $b} keys %$dataset )
    {
        my $id  = $dataset->{$ctr}->{'id'};
        my $seq = $dataset->{$ctr}->{'sequence'};

        if ( $id !~ /^PDBID/ )
        {
            if ( $seqctr < $SEQLIMIT )
            {
                print OUTFILE ">$id\n$seq\n";
                $seqctr++;
            }
        }
        else
        {
            print OUTFILE ">$id\n$seq\n";
        }
    }

    close OUTFILE;
}
else
{
    backticks("cp -f $TMP/$seekid.out $OUTPUTFILE");
}


cleanup();



####################
####################


sub parseFasta
{
    my $inputFile = shift;

	my %entries = ();
	return \%entries unless -e $inputFile;

	# read file
	open INFILE, "<$inputFile" or return \%entries;

	my $id = "";
	my $sequence = "";
	my $ctr = 0;

	while( my $line = <INFILE> )
	{
    	chomp $line;
    	$line =~ s/[\r\n]//g;

        if( $line =~ /^>(.+)/ )
    	{

    		if ( $id ne "" && $sequence ne "" )
    		{
    		    $entries{$ctr}{'id'} = $id;
    		    $entries{$ctr}{'sequence'} = $sequence;
    		    $ctr++;

    		}

    		$id = $1;
    		#$id =~ s/\W/_/g;   #safename
    		$sequence = "";
    	}
    	else
    	{
    		next if $id eq "";
    		$line =~ s/\s//g;
    		$sequence .= $line;
    	}
	}

	if ( $id ne "" && $sequence ne "" )
	{
	    $entries{$ctr}{'id'} = $id;
	    $entries{$ctr}{'sequence'} = $sequence;
	    $ctr++;
	}


	close INFILE;

	return \%entries;

}

sub parseResponse
{
    my $response = shift;
    my $status = "";
    my $seekid = "";

    if ( $response =~ /^([^\s:]+):([^\s:]+)$/ )
    {
        $seekid = $1;
        $status = $2;
    }

    return ($status, $seekid);
}


sub parseJobQueryResponse
{
    my $response = shift;
    my $jobs = 100;

    if ( $response =~ /^(\d+)$/ )
    {
        $jobs = $1;
    }

    return $jobs;
}



sub extractchecksum
{
    my $infile = shift;
    my %dataset = ();

    open CSUM, "tar -zxf $infile -O|" or return \%dataset;

    while(<CSUM>)
    {
        chomp;
        if ( /^(\S+)\s+(\S+)$/ )
        {
            $dataset{$2} = $1;
        }

    }

    close CSUM;

    return \%dataset;

}


sub parseError
{
    my $response = shift;

    #"error":"Invalid number of inputs found."
    my $errorstr = ( $response =~ /\"error\"\s*:\s*\"([^\"]+)\"/ ) ? $1 : $response;
    return $errorstr;
}


sub getchecksum
{
    my $infile = shift;

    # md5 binary check
    my $MD5BIN = "";

    if ( -x "/usr/bin/md5sum" )
    {
        $MD5BIN = "/usr/bin/md5sum";
    }
    elsif ( -x "/sbin/md5" )
    {
        $MD5BIN = "/sbin/md5 -q";
    }

    return "" if $MD5BIN eq "";


    my $checksum = "";
    open MD5EXE, "$MD5BIN $infile|" or return "";

    while(<MD5EXE>)
    {
        if (/^(\S+)\s+(\S+)$/)
        {
            $checksum = $1;
            last;
        }
        elsif (/^(\S+)$/)
        {
            $checksum = $1;
            last;
        }
    }

    close MD5EXE;

    return $checksum;

}


sub backticks
{
    my $command = shift;

    `$command`;
    return ($? == -1) ? 0 : 1;
}


sub bail
{
    my $str = shift;
    my $status = shift;

    #0 for success and 1 for error
    $status = 1 unless defined;

    print STDERR "$str\n" if defined $str;

    cleanup();

    exit($status);
}


sub cleanup
{
    return if ($TMP eq "" || !-d $TMP);

    opendir(MAINDIR, $TMP);
    my @files = readdir(MAINDIR);
    closedir(MAINDIR);

    foreach my $file (@files)
    {
        unlink "$TMP/$file" if -e "$TMP/$file";
    }

    remove_tree($TMP);

}


sub help
{
    my $str = shift;

    print <<'HELPME';

USAGE
 ./seekquencer_premafft.pl -inp <INFILE> -out <OUTFILE> [-str] [-seq] [-lim <count>] [-mod <mafftash|mafftash-split|fasta>] [-pre] [-run <thread|normal>] [-trd <count>]
 ./seekquencer_premafft.pl -idf <LISTFILE> -seqf <SEQFASTA> -out <OUTFILE> [-str] [-seq] [-lim <count>] [-mod <mafftash|mafftash-split|fasta>] [-pre] [-run <thread|normal>] [-trd <count>]


PARAMETERS
  -inp <INFILE>
     INFILE is a FASTA-formatted file.
     PDB entries are written as:
        >PDBID
        [5-character pdbid+chain]

     While sequence entries are written as:
        >[id]
        [sequence]

  -idf <LISTFILE>
     IDLISTFILE is a file containing a list of pdbids
     pdbids should be a 5-character pdbid + chain

  -seqf <SEQFASTA>
     SEQFASTA is a fasta file
     entries are written as:
        >[id]
        [sequence]

  -out <OUTFILE>
     Results are writen to a file named OUTFILE.

  -str
     Only structures will be collected by Seekquencer.
     If -str and -seq are not set, both structures and sequences will be collected by Seekquencer.

  -seq
     Only sequences will be collected by Seekquencer
     If -str and -seq are not set, both structures and sequences will be collected by Seekquencer.

  -lim <count>
     this sets the maximum number of sequence homologs collected. Default value: 100

  -pre
     When -str is set, this will compare all structures against all using pdp-ash.
     This would ensure that all structures collected are matching.
     All structures that do not match will be removed.

  -noin
     When set, inputs will not be included in the output

  -mod <mafftash|mafftash-split|fasta>
     Defines the output format.
     mafftash (default) will print a mafftash-formatted fasta file
     mafftash-split will make 2 files separating the structures (OUTFILE.str) from sequences (OUTFILE.seq)
     fasta will print a regular fasta file

  -run <thread|normal>
    thread will run simultaneous jobs during blast queries (faster but takes more nodes)
    normal will run sequential blast queries (slower but takes less nodes)
    Default value: normal

  -trd <count>
    if -run <thread> is defined, this sets the number of parallel jobs to run. Default value: 3


HELPME

    bail($str);
}



