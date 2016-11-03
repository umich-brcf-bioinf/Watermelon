#!/usr/bin/perl

### Parts of the code came from Perl tutorial on the web
### Usage: perl /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/send_email.pl <message> [starttime] 
use strict;
use warnings;
use POSIX;
use constant DATE => strftime("%a, %B %e, %Y", localtime);
use constant TIME => strftime("%I:%M %p", localtime);
use MIME::Lite;
my $message = $ARGV[0];
my $start;
my $error_log;
if ($ARGV[1]) {
	$start = $ARGV[1];
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

my $username = getpwuid( $< );
#print "User is $username\n";

my $to = "${username}\@umich.edu";
my $cc = "mpande\@umich.edu";
my $from = "mpande\@umich.edu";
my $subject = 'RNA-seq run';

my $msg = MIME::Lite->new(
                 From     => $from,
                 To       => $to,
                 Cc       => $cc,
                 Subject  => $subject,
                 Data     => $message
                 );
                 
$msg->send;
print "Email Sent Successfully\n";

