#!/usr/bin/perl

### Parts of the code came from Perl tutorial on the web
### Usage: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/send_email.pl <message> <cc> <error_log> [starttime]  
use strict;
use warnings;
use POSIX;

use constant DATE => strftime("%a, %B %e, %Y", localtime);
use constant TIME => strftime("%I:%M %p", localtime);
use MIME::Lite;
my $message = $ARGV[0];
my $collaborator = $ARGV[1];
my $error_log = $ARGV[2];
my $start;

if ($ARGV[2]) {
	$start = $ARGV[3];
} 
#else {$start = 1395765508;} ## on March 25, 2014 12:38 PM, seconds since epoch (00:00:00 UTC, January 1, 1970)
my $end = time();

if ($start) {
	my $runtime = $end - $start;
	my $days = int($runtime/(60*60*24));
	my $hours = ($runtime/(60*60))%24;
	my $minutes = ($runtime/60)%60;
	my $seconds = $runtime%60;
	#print "$days days, $hours hours, $minutes minutes and $seconds seconds.\n";
	$hours = $hours + ($days*24);
	$message = $message . DATE . " at " . TIME . ".\n" . "Total run time: $hours:$minutes:$seconds.\n";
}
if (-e $error_log) {
	$message = $message . "\nSee error log for error message/s. See run log for details.";
}
my $username = getpwuid( $< );
#print "User is $username\n";

my $to = "${username}\@umich.edu";
my $bcc = "mpande\@umich.edu";
my $from = "mpande\@umich.edu";
my $subject = 'RNA-seq run';
my $cc = "";
if ($collaborator && $collaborator ne "NA") {
	my @collaborators = split(",", $collaborator);
	my $c1 = shift(@collaborators);
	if ($c1 !~ /@/) {
		$c1 = $c1 . "\@umich.edu";
	}
	$cc = $c1;
	foreach my $c (@collaborators) {
		if ($c !~ /@/) {
			$c = $c . "\@umich.edu";
		}
		$cc = $cc . ", $c";
	}
}
#print "\$cc is $cc\n";
#=pod
my $msg = MIME::Lite->new(
                 From     => $from,
                 To       => $to,
                 Cc       => $cc,
		   Bcc       => $bcc,	
                 Subject  => $subject,
                 Data     => $message
                 );
                 
$msg->send;
print "Email Sent Successfully\n";
#=cut
