#!/usr/bin/env perl

#print $ARGV[0]." ".$ARGV[1]." ".$#ARGV."\n";



if(($#ARGV != 1))
  {
    print "ERROR: usage: perl applyRAxML2AllFilesInDirectory.pl directoryName \"raxmlCommandLine\" \n";
    exit;
  }

$directory = $ARGV[0];
$commandLine = $ARGV[1];
@fileList;

opendir (DIR, $directory) or die $!;

while (my $file = readdir(DIR)) 
{
    #print "$file\n";
    push(@fileList, $file); 
}

chdir DIR;

closedir(DIR);



foreach (@fileList)
{
    $string = $_;
   
    if($string !~ /^[.]/ && $string !~ /.+~/)
    {
	$cmd = $commandLine." -n Analysis.".$string." -s ".$string."\n";
	print "$cmd";
	system($cmd);
    }

   
}


