#!/usr/local/bin/perl

#print $ARGV[0]." ".$ARGV[1]." ".$#ARGV."\n";

#adapt as required, modify to read ./raxmHPC if raxml is not in your path



$raxmlExecutable = "raxmlHPC-AVX";
# $raxmlExecutable = "raxmlHPC-SSE3";
# $raxmlExecutable = "raxmlHPC-PTHREADS -T 4";
# $raxmlExecutable = "raxmlHPC-PTHREADS-AVX -T 4";
# $raxmlExecutable = "raxmlHPC-PTHREADS-SSE3 -T 4";

if(($#ARGV < 0) || ($#ARGV > 1))
  {
    print "ERROR: usage: \"perl ProteinModelSelection.pl alignmentFileName [modelFileName] \"\n";
    exit;
  }

$alignmentName = $ARGV[0];

$UNLIKELY = -1.0E300;

sub getLH
  {
    my $fileID = $_[0];  
    open(CPF, $fileID);
    my @lines = <CPF>;	   
    my $numIT = @lines;   	
    my $lastLH = pop(@lines);  
    my $k = index($lastLH, '-');   
    my $LH = substr($lastLH, $k);     
    close(CPF);	
    return $LH;
  }

sub getTIME
  {
    my $fileID = $_[0];  
    open(CPF, $fileID);
    my @lines = <CPF>;	   	
    my $numIT = @lines;   	
    my $lastLH = pop(@lines);  
    my $k = index($lastLH, '-');   
    my $TIME = substr($lastLH, 0, $k-1);    
    close(CPF);  
    return $TIME;
  }


@AA_Models = ("DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", 
	      "LG", "MTART", "MTZOA", "PMB", "HIVB", "HIVW", "JTTDCMUT", "FLU",
	      "DAYHOFFF", "DCMUTF", "JTTF", "MTREVF", "WAGF", "RTREVF", "CPREVF", "VTF", "BLOSUM62F", 
	      "MTMAMF", "LGF", "MTARTF", "MTZOAF", "PMBF", "HIVBF", "HIVWF", "JTTDCMUTF", "FLUF");



if($#ARGV == 1)
  {
    print "Splitting up multi-gene alignment\n";
    $partition =  $ARGV[1];
    $cmd = $raxmlExecutable." -f s -m PROTCATJTT -O -p 12345 -s ".$alignmentName." -q ".$partition." -n SPLIT_".$alignmentName." \> SPLIT_".$alignmentName."_out";
    system($cmd);
    $count = 0;
    while(open(CPF, $alignmentName.".GENE.".$count))
      {
	close CPF;
	print "PARTITION: ".$count."\n";
	#print "perl ProteinModelSelection.pl ".$alignmentName.".GENE.".$count;
	system("perl ProteinModelSelection.pl ".$alignmentName.".GENE.".$count);
	$count = $count + 1;
      }
  }
else
  {
    #print "Determining AA model data\n";
    #print "Computing randomized stepwise addition starting tree number :".$i."\n";
    $cmd = $raxmlExecutable." -O -y -p 12345 -m PROTCATJTT -s ".$alignmentName." -n ST_".$alignmentName." \> ST_".$alignmentName."_out";
    system($cmd);
    
    $numberOfModels = @AA_Models;
    
    for($i = 0; $i < $numberOfModels; $i++)
      {
	$aa = "PROTGAMMA".$AA_Models[$i];
	$cmd = $raxmlExecutable." -O -f e -m ".$aa." -s ".$alignmentName." -t RAxML_parsimonyTree.ST_".$alignmentName." -n ".$AA_Models[$i]."_".$alignmentName."_EVAL \> ".$AA_Models[$i]."_".$alignmentName."_EVAL.out\n";  
	#print($cmd);
	system($cmd);
      }
    
   
    for($i = 0; $i < $numberOfModels; $i++)
      {
	$logFileName = "RAxML_log.".$AA_Models[$i]."_".$alignmentName."_EVAL";
	#print $logFileName."\n";
	$lh[$i] = getLH($logFileName);
      }
    
    $bestLH = $UNLIKELY;
    $bestI = -1;
    
    for($i = 0; $i < $numberOfModels; $i++)
      {
	#print "Model: ".$AA_Models[$i]." LH: ". $lh[$i]."\n";
	if($lh[$i] > $bestLH)
	  {
	    $bestLH = $lh[$i];
	    $bestI = $i;
	  }
      }
    
    print "Best Model : ".$AA_Models[$bestI]."\n\n";
    $bestModel = $AA_Models[$bestI]; 
  }
    
