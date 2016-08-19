#!/usr/bin/perl
#use String::Scanf; 
open FILE, "<", "id_squid" or die $!;

while (<FILE>)
{
        if($_ =~ m/(F.+P.+L.+S)(.+)(=.+)\n/)
	{	

		print $1;
#		$num = sscanf("%f",$2);
		if($2>0.1)
		{
			print " $2 ";
		}
		else
		{
			print " 0.100 ";
		}
		print $3;
		print "\n";
	}
	else
	{
		print $_;
	}
}


